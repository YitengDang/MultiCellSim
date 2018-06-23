function cells = init_random_cells_montecarlo(N, p)
    % Monte Carlo algorithm to generate random configuration with fixed N, p
    % Do not simulate if p=0 or p=1
    if p==0
        cells = zeros(N,1);
        return
    elseif p==1
        cells = ones(N,1);
        return
    end
    % If p is not 0 or 1
    mcsteps = 10^3;
    cells = p*ones(N, 1);
    r = randi(N, mcsteps, 2);
    for i=1:mcsteps
        delta = (2*rand()-1)/2; % -0.5 <= delta <= 0.5
        cell1 = cells(r(i, 1)) + delta;
        cell2 = cells(r(i, 2)) - delta;
        if abs(2*cell1-1)<=1 && abs(2*cell2-1)<=1 % if both cells still in range [0, 1]
            cells(r(i, 1)) = cell1;
            cells(r(i, 2)) = cell2;
        end
    end
end