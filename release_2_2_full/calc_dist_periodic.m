function dist = calc_dist_periodic(x, y, Lx, Ly)
    % calculates the distances between cells for periodic boundaries
    %x = pos(:, 1);
    %y = pos(:, 2);
    [x1,x2] = meshgrid(x, x);
    [y1,y2] = meshgrid(y, y);
    dx = mod(abs(x1 - x2), Lx);
    dx = min(dx, Lx-dx);
    dy = mod(abs(y1 - y2), Ly);
    dy = min(dy, Ly-dy);
    dist = sqrt(dx.^2 + dy.^2);
end
