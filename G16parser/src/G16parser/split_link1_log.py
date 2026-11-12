import os
import re

from ._common import read_lines

_CHK_RE = re.compile(r"^\s*%chk\s*=\s*(\S+)", re.IGNORECASE)
_DEFAULT_PATTERN = r"scan along X: N = (-?\d+)"


def g16_split_link1_log(logfile, naming="auto", pattern=_DEFAULT_PATTERN,
                         prefix="", outdir=""):
    """Splits a Gaussian log of chained (--Link1--) jobs into one log per job.

    The Gaussian output format for --Link1-- chaining does not differ
    between Gaussian 09 and Gaussian 16, so this single function covers
    both (matching MATLAB's split_link1_log.m, shared identically between
    the G09/ and G16/ folders).

        g16_split_link1_log('NMA_WB97XD_DEF2TZVP_TD_ALL.LOG')
        g16_split_link1_log('NMA_wB97XD_def2TZVP_scan_all.log')        # field scans
        files = g16_split_link1_log(..., naming='chk')

    Two ways of naming the output logs (`naming`):
        'auto'  (default) field-scan title if found in the first job
                ("... scan along X: N = ...", as written by make_field_scan.m),
                otherwise the %chk name of each job.
        'field' <prefix>_Fp020.log, <prefix>_Fm040.log, ... from the title
                (old behaviour; raises if the title is missing).
        'chk'   <name of the job's .chk file>.log, e.g.
                NMA_wB97XD_def2TZVP_TD_trans.log, NMA_wB97XD_IP_EA ... .log

    The internal Link1 step that Gaussian creates for the Freq part of an
    Opt+Freq job carries the same title / %chk (or none) and is kept
    together with its optimization.

    Parameters
    ----------
    naming  : 'auto' | 'field' | 'chk'          (default 'auto')
    pattern : regex with one capturing group for the field value N in the
              title (default 'scan along X: N = (-?\\d+)')
    prefix  : field naming only; default = log name without a trailing
              '_scan_all'
    outdir  : default = folder of the input log

    Returns
    -------
    files : list[str] with the full paths of the logs written.
    """
    if naming not in ("auto", "field", "chk"):
        raise ValueError(
            f"g16_split_link1_log: naming must be 'auto', 'field' or 'chk', got {naming!r}."
        )

    d, b = os.path.split(logfile)
    base, _ext = os.path.splitext(b)
    if not outdir:
        outdir = d if d else os.getcwd()
    if not prefix:
        prefix = re.sub(r"_scan_all$", "", base, flags=re.IGNORECASE)
    os.makedirs(outdir, exist_ok=True)

    lines = read_lines(logfile)
    nlines = len(lines)

    # ---------- segments between "Entering Link 1" lines ----------
    starts = [i for i, ln in enumerate(lines) if "Entering Link 1" in ln]
    if not starts:
        raise ValueError(f"No \"Entering Link 1\" line found in {logfile}.")
    starts[0] = 0   # keep the file header with the first job
    ends = [s - 1 for s in starts[1:]] + [nlines - 1]
    nseg = len(starts)

    # ---------- key of each segment: field value N or %chk name ----------
    field_re = re.compile(pattern)
    seg_n = [None] * nseg
    seg_chk = [""] * nseg
    for s in range(nseg):
        blk = lines[starts[s]:ends[s] + 1]
        for ln in blk:
            m = field_re.search(ln)
            if m:
                seg_n[s] = float(m.group(1))
                break
        for j, ln in enumerate(blk):
            m = _CHK_RE.match(ln)
            if m:
                # Gaussian wraps long Link0 lines at 80 columns: join the
                # continuation lines until the name is complete (ends in .chk)
                cname = m.group(1).strip()
                k = j + 1
                while (not cname.lower().endswith(".chk")
                       and k < len(blk) and k <= j + 4):
                    nxt = blk[k].strip()
                    if nxt == "" or nxt.startswith(("-", "%", "#")):
                        break
                    cname += nxt
                    k += 1
                cname = cname.replace("\\", "/")
                seg_chk[s] = os.path.splitext(os.path.basename(cname))[0]
                break

    resolved = naming
    if resolved == "auto":
        resolved = "field" if seg_n[0] is not None else "chk"

    if resolved == "field":
        if seg_n[0] is None:
            raise ValueError(
                "Title pattern not found in the first job. Check 'pattern', "
                "or use naming='chk'."
            )
        for s in range(1, nseg):
            if seg_n[s] is None:
                seg_n[s] = seg_n[s - 1]
        key = seg_n
    else:
        if seg_chk[0] == "":
            raise ValueError(
                "No chk line (%chk=...) found in the first job: use naming='field' "
                "or add a %chk line to the input."
            )
        for s in range(1, nseg):
            if seg_chk[s] == "":
                seg_chk[s] = seg_chk[s - 1]
        key = seg_chk

    # ---------- group consecutive segments with the same key ----------
    g_start = [0] + [i for i in range(1, nseg) if key[i] != key[i - 1]]
    g_end = [i - 1 for i in g_start[1:]] + [nseg - 1]

    print(f"{'File':<46} {'Job':>5} {'Normal':>8} Status")
    written = []
    files = []
    input_log_abspath = os.path.abspath(os.path.join(outdir, b))
    for g, (gs, ge) in enumerate(zip(g_start, g_end), start=1):
        a, z = starts[gs], ends[ge]
        if resolved == "field":
            n = key[gs]
            name = (f"{prefix}_Fp{int(n):03d}" if n >= 0
                     else f"{prefix}_Fm{int(-n):03d}")
        else:
            name = str(key[gs])

        base_name, k = name, 1
        while name in written:
            k += 1
            name = f"{base_name}_{k}"
        written.append(name)

        out_path = os.path.join(outdir, name + ".log")
        if os.path.abspath(out_path) == input_log_abspath:
            name = name + "_job"             # never overwrite the input log
            out_path = os.path.join(outdir, name + ".log")
        files.append(out_path)

        blk = lines[a:z + 1]
        with open(out_path, "w", encoding="utf-8") as f:
            f.write("\n".join(blk) + "\n")

        n_norm = sum("Normal termination" in ln for ln in blk)
        n_link = ge - gs + 1                 # Gaussian link steps in this job
        if any("Error termination" in ln for ln in blk):
            status = "ERROR termination"
        elif n_norm >= n_link:
            status = "OK"
        elif n_norm > 0:
            status = f"{n_norm} of {n_link} steps terminated"
        else:
            status = "incomplete / still running"
        if any("Stationary point found" in ln for ln in blk):
            status += ", stationary point"
        print(f"{name + '.log':<46} {g:>5} {n_norm:>8} {status}")

    print(f"{len(g_start)} job logs written to {outdir}")
    return files
