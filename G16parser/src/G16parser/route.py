import re

from ._common import read_lines

_HASH_RE = re.compile(r"^#")


def g16_route(filename, lines=None):
    """Extracts the route section from a Gaussian 16 .out/.log file.

    Collects the lines between the two '----' separators that follow the
    first '#' line (the standard Gaussian route block).

    Parameters
    ----------
    filename : str
    lines : list[str], optional
        Pre-read file lines, to skip re-reading the file when it has
        already been read elsewhere (see g16_read_all).

    Returns
    -------
    route : str — full route section string, on a single line
    """
    if lines is None:
        lines = read_lines(filename)

    route_lines = []
    in_route = False
    found_first_sep = False

    for raw_ln in lines:
        ln = raw_ln.strip()

        is_sep = bool(ln) and set(ln) == {"-"} and len(ln) >= 20

        if is_sep:
            if in_route:
                break
            found_first_sep = True
            continue

        if found_first_sep and not in_route:
            if _HASH_RE.match(ln):
                in_route = True
                # Keep only the leading indent stripped; preserve any
                # trailing whitespace exactly as printed (see join
                # comment below).
                route_lines.append(raw_ln.lstrip())
            else:
                found_first_sep = False
            continue

        if in_route:
            route_lines.append(raw_ln.lstrip())

    if not route_lines:
        raise ValueError(f"g16_route: route section not found in {filename}")

    # Gaussian wraps the route echo at a fixed column with no regard for
    # word boundaries, so a keyword can be split mid-word across two
    # lines (e.g. "nosym" / "m cphf=..." for "nosymm cphf=..."). Joining
    # with an inserted space (the previous behaviour) corrupted every
    # such keyword. Gaussian never leaves a genuine trailing space before
    # the forced wrap, so concatenating with no separator reconstructs
    # the original text correctly whether the wrap fell mid-word or
    # between two words (any real separating space is part of the line
    # content itself, not the join).
    route = re.sub(r"\s+", " ", "".join(route_lines).strip())
    print(f"Route: {route}")
    return route
