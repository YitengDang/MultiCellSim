function dist = calc_dist(x, y, Lx, Ly, bc_periodic)
    % calculates the distances between cells for periodic boundaries
    %x = pos(:, 1);
    %y = pos(:, 2);
    % bc_periodic = 2x1 array with values in {0,1}
    % boundary conditions -> (0,0) = fixed, (1,1) = periodic, (1,0)
    % = periodic in x-axis only, (0, 1) = periodic in y-axis only
    
    [x1,x2] = meshgrid(x, x);
    [y1,y2] = meshgrid(y, y);
    dx = mod(abs(x1 - x2), Lx);
    dy = mod(abs(y1 - y2), Ly);
    if bc_periodic(1)
        dx = min(dx, Lx-dx);
    end
    if bc_periodic(2)
        dy = min(dy, Ly-dy);
    end
    dist = sqrt(dx.^2 + dy.^2);
end
