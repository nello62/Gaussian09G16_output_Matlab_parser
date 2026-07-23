import re

import numpy as np

from ._common import Struct, read_lines

_SEPARATOR_RE = re.compile(r"^-{4,}$")
_LINK0_RE = re.compile(r"^%(\w+)\s*=\s*(.+)$")
_CHARGE_MULT_RE = re.compile(r"^(-?\d+)\s+(\d+)")
_ATOM_SYMBOL_RE = re.compile(r"^([A-Za-z]{1,2})")


def g16_read_input(filename):
    """Reads a Gaussian input file (.gjf/.com/.in) and extracts link0
    commands, route section, title, charge/multiplicity, and the starting
    Cartesian geometry.

    The Gaussian input file format does not differ between Gaussian 09 and
    Gaussian 16, so this single function covers both (matching MATLAB's
    G_read_input.m, shared identically between the G09/ and G16/ folders).

    Returns
    -------
    ginp : Struct with fields
        symbols   list[str]           atomic symbols
        xyz       np.ndarray (N, 3)   starting Cartesian coordinates (Angstrom)
        Natoms    int
        charge    int
        mult      int
        title     str      title/comment line(s), joined with a space
        route     str      full route section, single line
        method    str      method parsed from the route (e.g. 'b3lyp'), '' if not found
        basis     str      basis set parsed from the route (e.g. '6-311+g(d,p)'), '' if not found
        chk       str      %chk link0 value, '' if absent
        mem       str      %mem link0 value, '' if absent
        nproc     str      %nprocshared/%nproc link0 value, '' if absent
        filename  str      source file path

    ginp has the same .symbols/.xyz/.Natoms/.filename fields as the Struct
    returned by g16_structure, so it can be passed directly to
    g16_draw_molecule or g16_get_bond_length without any conversion.

    Limitation: only Cartesian-coordinate geometry blocks are supported
    (not Z-matrix input).
    """
    lines = read_lines(filename)

    # Drop pure comment lines ('!' as the first non-blank character), which
    # Gaussian allows anywhere in an input file and which carry no structural
    # meaning for this parser.
    lines = [ln for ln in lines if not ln.strip().startswith("!")]
    n = len(lines)

    def is_blank(idx):
        return idx < n and lines[idx].strip() == ""

    def is_separator(idx):
        return idx < n and bool(_SEPARATOR_RE.match(lines[idx].strip()))

    # -------------------------------------------------------------------
    # Link0 (%) commands
    # -------------------------------------------------------------------
    chk = mem = nproc = ""
    i = 0
    while i < n:
        ln = lines[i].strip()
        if ln == "":
            i += 1
            continue
        if not ln.startswith("%"):
            break
        m = _LINK0_RE.match(ln)
        if m:
            key = m.group(1).lower()
            val = m.group(2).strip()
            if key == "chk":
                chk = val
            elif key == "mem":
                mem = val
            elif key in ("nprocshared", "nproc"):
                nproc = val
        i += 1

    # Skip a decorative separator line before the route, if present
    while is_separator(i):
        i += 1
    while is_blank(i):
        i += 1

    # -------------------------------------------------------------------
    # Route section: from the first '#' line to the next blank line
    # -------------------------------------------------------------------
    route_lines = []
    first_trimmed = lines[i].strip() if i < n else ""
    if first_trimmed.startswith("#"):
        while i < n:
            tln = lines[i].strip()
            if tln == "" or is_separator(i):
                # A blank line or a closing separator both end the route:
                # some generated inputs (e.g. G09_restart) bracket the route
                # between two separator lines with no blank line before the
                # title that follows, so a separator must terminate the
                # route immediately rather than being skipped over.
                i += 1
                break
            route_lines.append(tln)
            i += 1
    route = " ".join(route_lines).strip()

    method = basis = ""
    rtoks = route.split()
    mb_idx = next((k for k, t in enumerate(rtoks) if "/" in t), None)
    if mb_idx is not None:
        parts = rtoks[mb_idx].split("/")
        if len(parts) >= 2:
            method, basis = parts[0], parts[1]

    # Skip a decorative separator line after the route, and blank lines
    while i < n and (lines[i].strip() == "" or is_separator(i)):
        i += 1

    # -------------------------------------------------------------------
    # Title (one or more non-blank lines, until the next blank line)
    # -------------------------------------------------------------------
    title_lines = []
    while i < n and lines[i].strip() != "":
        title_lines.append(lines[i].strip())
        i += 1
    titlestr = " ".join(title_lines)
    i += 1   # skip the blank line after the title

    # -------------------------------------------------------------------
    # Charge / multiplicity line
    # -------------------------------------------------------------------
    charge = mult = None
    while is_blank(i):
        i += 1
    if i < n:
        m = _CHARGE_MULT_RE.match(lines[i].strip())
        if m:
            charge = int(m.group(1))
            mult = int(m.group(2))
        i += 1

    # -------------------------------------------------------------------
    # Geometry: "Symbol  X  Y  Z" Cartesian lines, until a blank line or EOF
    # -------------------------------------------------------------------
    symbols = []
    xyz_rows = []
    while i < n:
        ln = lines[i].strip()
        if ln == "":
            break
        parts = ln.split()
        if len(parts) >= 4:
            try:
                x, y, z = float(parts[-3]), float(parts[-2]), float(parts[-1])
            except ValueError:
                x = y = z = None
            if x is not None and not any(np.isnan([x, y, z])):
                symtok = _ATOM_SYMBOL_RE.match(parts[0])
                if symtok:
                    sym = symtok.group(1)
                    sym = sym[0].upper() + sym[1:].lower() if len(sym) > 1 else sym.upper()
                    symbols.append(sym)
                    xyz_rows.append((x, y, z))
        i += 1

    if not symbols:
        raise ValueError(
            f"g16_read_input: no Cartesian geometry block found in {filename} "
            "(Z-matrix input is not supported)."
        )

    print(f"\n── g16_read_input: {filename} ──")
    print(f"  Route  : {route}")
    if method:
        print(f"  Method/Basis : {method} / {basis}")
    print(f"  Charge = {charge}   Multiplicity = {mult}")
    print(f"  {len(symbols)} atoms\n")

    return Struct(
        symbols=symbols,
        xyz=np.array(xyz_rows, dtype=float),
        Natoms=len(symbols),
        charge=charge,
        mult=mult,
        title=titlestr,
        route=route,
        method=method,
        basis=basis,
        chk=chk,
        mem=mem,
        nproc=nproc,
        filename=filename,
    )
