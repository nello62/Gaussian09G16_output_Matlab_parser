import os
import re
import warnings

from ._common import Struct


def g16_gaussian_version(filename):
    """Detects which Gaussian version/revision produced a .out/.log/.fchk file.

    For .out/.log files this reads the "Gaussian NN, Revision X.YY," citation
    line that Gaussian prints near the top of every output file (works for
    any NN: 09, 16, ...).

    .fchk files do NOT store this information at all — the formatted
    checkpoint format has no version/provenance field. In that case this
    looks for a sibling .log/.out file with the same base name in the same
    folder and reads the version from there instead. If no sibling file is
    found, .major/.revision/.full come back None and a warning is issued.

    Returns a Struct with fields:
        major     int    Gaussian major version, e.g. 9, 16 (None if unknown)
        revision  str    revision string, e.g. 'A.02', 'C.01' ('' if unknown)
        full      str    citation line as printed (e.g. 'Gaussian 16, Revision C.01')
        source    str    'out/log' | 'fchk-sibling:<path>' | 'unknown'
        filename  str
    """
    if not os.path.isfile(filename):
        raise FileNotFoundError(f"g16_gaussian_version: file not found: {filename}")

    ext = os.path.splitext(filename)[1].lower()

    if ext == ".fchk":
        major, revision, full = _read_version(filename)
        if major is not None:
            source = "out/log"  # unexpected for a genuine .fchk, kept for completeness
        else:
            sibling = _find_sibling(filename)
            if sibling:
                major, revision, full = _read_version(sibling)
                source = f"fchk-sibling:{sibling}" if major is not None else "unknown"
            else:
                source = "unknown"
    else:
        major, revision, full = _read_version(filename)
        source = "out/log" if major is not None else "unknown"

    if major is None:
        warnings.warn(
            f'Could not determine the Gaussian version for {filename} '
            f'(no "Gaussian NN, Revision ..." line found, and no readable '
            f'sibling .log/.out file).'
        )

    gv = Struct(major=major, revision=revision, full=full, source=source, filename=filename)

    if major is None:
        print(f"{filename} -> unknown Gaussian version")
    else:
        tail = f" [from {source}]" if source != "out/log" else ""
        print(f"{filename} -> {full} (revision {revision}){tail}")

    return gv


_VERSION_RE = re.compile(r"Gaussian\s+(\d+)\s*,\s*Revision\s+([\w.]+)")


def _read_version(fpath):
    with open(fpath, "r", encoding="latin-1", errors="replace") as f:
        for n, line in enumerate(f):
            m = _VERSION_RE.search(line)
            if m:
                major = int(m.group(1))
                revision = m.group(2)
                full = re.sub(r",\s*$", "", line.strip())
                return major, revision, full
            if n > 500:
                break
    return None, "", ""


def _find_sibling(fchk_path):
    folder, base = os.path.split(fchk_path)
    base = os.path.splitext(base)[0]
    for ext in (".log", ".out", ".LOG", ".OUT"):
        cand = os.path.join(folder, base + ext)
        if os.path.isfile(cand):
            return cand
    return ""
