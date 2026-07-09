import re
import warnings

import numpy as np

from ._common import Struct, read_lines, z_to_symbol


def g16_structure(filename, orientation="auto", step="last", lines=None):
    """Extracts the molecular geometry from a Gaussian 16 .out/.log file.

    Parameters
    ----------
    filename : str
    orientation : {'auto', 'standard', 'input'}
        'auto' (default) uses Standard orientation if present, otherwise
        falls back to Input orientation.
    step : 'last' (default) | 'first' | int (1-based)
    lines : list[str], optional
        Pre-read file lines, to skip re-reading the file when it has
        already been read elsewhere (see g16_read_all).

    Returns
    -------
    mol : Struct with fields
        symbols   list[str]           atomic symbols
        xyz       np.ndarray (N, 3)   Cartesian coordinates (Angstrom)
        Z         np.ndarray (N,)     atomic numbers
        Natoms    int
        step      int                 index of the extracted step (1-based)
        n_steps   int                 total geometry blocks in the file
        orientation str               'Standard orientation' | 'Input orientation'
        filename  str
    """
    orientation = orientation.lower()
    if orientation not in ("auto", "standard", "input"):
        raise ValueError("g16_structure: orientation must be 'auto', 'standard', or 'input'.")

    if lines is None:
        lines = read_lines(filename)

    if orientation == "auto":
        has_std = any(re.search(r"Standard orientation\s*:", ln, re.IGNORECASE) for ln in lines)
        ori_label = "Standard orientation" if has_std else "Input orientation"
    elif orientation == "standard":
        ori_label = "Standard orientation"
    else:
        ori_label = "Input orientation"

    pat = re.compile(ori_label + r"\s*:", re.IGNORECASE)
    block_starts = [i for i, ln in enumerate(lines) if pat.search(ln)]

    if not block_starts and orientation == "standard":
        warnings.warn("g16_structure: Standard orientation not found, falling back to Input orientation.")
        ori_label = "Input orientation"
        pat = re.compile(ori_label + r"\s*:", re.IGNORECASE)
        block_starts = [i for i, ln in enumerate(lines) if pat.search(ln)]

    if not block_starts:
        raise ValueError(f'g16_structure: no "{ori_label}" block found in {filename}')

    n_blocks = len(block_starts)
    if isinstance(step, str):
        if step.lower() == "last":
            step_idx = n_blocks
        elif step.lower() == "first":
            step_idx = 1
        else:
            raise ValueError("g16_structure: step must be 'first', 'last', or an integer.")
    else:
        step_idx = round(step)
        if step_idx < 1 or step_idx > n_blocks:
            raise ValueError(f"g16_structure: step {step_idx} out of range [1, {n_blocks}].")

    header_line = block_starts[step_idx - 1]   # 0-based index of the header line
    data_start = header_line + 5                # skip sep, 2 column-header lines, sep

    symbols, xyz_rows, z_rows = [], [], []
    row_re = re.compile(r"^\s*(\d+)\s+(\d+)\s+(\d+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)")

    for ln in lines[data_start:]:
        s = ln.strip()
        if not s:
            continue
        if set(s) == {"-"}:
            break
        m = row_re.match(ln)
        if not m:
            continue
        znum = int(m.group(2))
        x, y, z = float(m.group(4)), float(m.group(5)), float(m.group(6))
        symbols.append(z_to_symbol(znum))
        xyz_rows.append((x, y, z))
        z_rows.append(znum)

    if not xyz_rows:
        raise ValueError(f"g16_structure: no atoms read from step {step_idx}.")

    mol = Struct(
        symbols=symbols,
        xyz=np.array(xyz_rows, dtype=float),
        Z=np.array(z_rows, dtype=int),
        Natoms=len(xyz_rows),
        step=step_idx,
        n_steps=n_blocks,
        orientation=ori_label,
        filename=filename,
    )
    return mol
