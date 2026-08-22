"""ramanbrowser -- click-to-select Raman/IR stick-spectrum browser.

Python port of G_Utility/RamanBrowser/G_raman_browser.m: plots intensity
vs. frequency as a stick ("stem") plot; clicking a stick selects the
nearest mode (highlighted), and a button renders that mode's 3D
displacement structure via g16_draw_mode.

Standalone G_Utility companion, not part of the G16parser package (same
status as MOSurface/python/'s own mosurface package) -- but it DOES
depend on G16parser for g16_structure/g16_nmodes/g16_draw_mode, the one
real dependency, same relationship MOSurface/python/ has.

    from ramanbrowser import g16_raman_browser
    g16_raman_browser('4-NTP.out')

Scope note (a real difference from the MATLAB original): G16parser only
targets Gaussian 16 output (documented, deliberate scope -- see the
toolbox manual's own note on this), so this port is G16-only, unlike
G_raman_browser.m, which auto-detects and also supports Gaussian 09 via
G09_structure/G09_nmodes/G09_draw_mode. There is no G09 equivalent
anywhere in the Python side of this project to dispatch to.
"""

from .raman_browser import g16_raman_browser

__all__ = ["g16_raman_browser"]
