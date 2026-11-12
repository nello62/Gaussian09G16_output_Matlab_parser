function G_relight(ax, v)
%G_RELIGHT  Re-position the lights of a 3D axes so that they follow the
%   current view (or a new one), e.g. after changing 'View' on a figure made
%   by G_draw_mo_surface, G_draw_density_surface, G_draw_esp_surface or
%   G_draw_cube_surface.
%
%   G_relight                 % current axes, current view
%   G_relight(ax)             % given axes (or an array of axes)
%   G_relight(ax, [az el])    % set the view first, then relight
%   G_relight([], [0 90])     % current axes, view from the top
%
%   The lighting scheme is the one used by the G_draw_*_surface functions:
%   a headlight (light at the camera position) plus a second light at
%   azimuth 45 deg and elevation 30 deg relative to the camera, so the
%   surfaces keep their shading gradient whatever the viewing direction.
%   Existing lights are re-used (the first two created in the axes);
%   if the axes has no lights, they are created with gouraud lighting.
%
%   Typical use:
%     mo = G16_fchk_read('NMA_wB97XD_augTZ_TD_trans.fchk');
%     G_draw_mo_surface(mo, 'HOMO');
%     G_relight(gca, [-60 40]);
%     exportgraphics(gca, 'HOMO.pdf', 'ContentType', 'vector');

if nargin < 1 || isempty(ax), ax = gca; end
for k = 1:numel(ax)
    a = ax(k);
    if ~isgraphics(a, 'axes')
        error('G_relight: input must be an axes handle.');
    end
    if nargin >= 2 && ~isempty(v)
        view(a, v);
    end
    L = findobj(a, 'Type', 'light');     % newest first
    L = flipud(L(:));                    % oldest (created first) first
    if isempty(L)
        lighting(a, 'gouraud');
        camlight(a, 'headlight');
        camlight(a, 45, 30);
    else
        camlight(L(1), 'headlight');
        if numel(L) >= 2
            camlight(L(2), 45, 30);
        else
            camlight(a, 45, 30);
        end
    end
end
end
