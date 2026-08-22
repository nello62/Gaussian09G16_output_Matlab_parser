"""MOSurface -- real-space MO/density/ESP isosurfaces and cube-file
load-and-replot, evaluated directly from Gaussian .fchk basis-set data.

Python port of the MATLAB G_Utility/MOSurface/ collection; see README.md
for the naming philosophy, dependency on G16parser, validation performed,
and known performance characteristics.

    from mosurface import g16_draw_mo_surface, g16_draw_density_surface
    import G16parser as g16
    data = g16.g16_fchk_read('molecule.fchk')
    g16_draw_mo_surface(data, 'HOMO')
"""

from .draw_mo_surface import g16_draw_mo_surface
from .draw_density_surface import g16_draw_density_surface
from .draw_esp_surface import g16_draw_esp_surface
from .draw_cube_surface import g16_draw_cube_surface

__all__ = [
    "g16_draw_mo_surface",
    "g16_draw_density_surface",
    "g16_draw_esp_surface",
    "g16_draw_cube_surface",
]
