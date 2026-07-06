function h = draw_rectangle(x1, y1, x2, y2, color, alpha)
% DRAW_RECTANGLE  Disegna un rettangolo riempito in un grafico X-Y.
%
%   h = draw_rectangle(x1, y1, x2, y2, color, alpha)
%
%   Input:
%     x1, y1  - coordinate del primo vertice (angolo in basso a sinistra)
%     x2, y2  - coordinate del vertice opposto (angolo in alto a destra)
%     color   - colore del rettangolo (stringa es. 'r', o tripla RGB es. [0.2 0.5 1])
%     alpha   - trasparenza: 0 = completamente trasparente, 1 = opaco
%
%   Output:
%     h       - handle all'oggetto patch creato
%
%   Esempio:
%     figure; hold on; axis([-1 5 -1 5]);
%     draw_rectangle(0, 0, 4, 3, [0.2 0.6 1], 0.4);
%     draw_rectangle(1, 1, 3, 2, 'r', 0.7);

    % Vertici del rettangolo (senso antiorario)
    xv = [x1, x2, x2, x1];
    yv = [y1, y1, y2, y2];

    % Disegna il rettangolo come patch
    h = patch(xv, yv, color, ...
              'FaceAlpha',  alpha, ...
              'EdgeColor',  color, ...
              'LineWidth',  1.5);
end