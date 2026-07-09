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
                route_lines.append(ln)
            else:
                found_first_sep = False
            continue

        if in_route:
            route_lines.append(ln)

    if not route_lines:
        raise ValueError(f"g16_route: route section not found in {filename}")

    route = " ".join(route_lines).strip()
    print(f"Route: {route}")
    return route
