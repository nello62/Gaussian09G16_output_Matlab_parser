import re

from ._common import read_lines

_RE = re.compile(r"Charge\s*=\s*([-\d]+)\s+Multiplicity\s*=\s*(\d+)")


def g16_charge_mult(filename, lines=None):
    """Extracts molecular charge and spin multiplicity from a Gaussian 16
    .out/.log file.

    Parameters
    ----------
    filename : str
    lines : list[str], optional
        Pre-read file lines (see g16_read_lines / g16_read_all), to skip
        re-reading the file when it has already been read elsewhere.

    Returns
    -------
    (charge, mult) : tuple[int, int]
    """
    if lines is None:
        lines = read_lines(filename)

    for line in lines:
        m = _RE.search(line)
        if m:
            charge = int(m.group(1))
            mult = int(m.group(2))
            print(f"Charge = {charge:+d}   Multiplicity = {mult}")
            return charge, mult

    raise ValueError(
        f'g16_charge_mult: "Charge = ... Multiplicity = ..." line not found in {filename}'
    )
