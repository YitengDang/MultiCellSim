function [cells_out, changed, mom, I] = ...
    update_cells_noise(cells, dist, Son, K, a0, Rcell, noise)
% Update cells using noise in a positive feedback loop with infinite hill
% coefficient

% Account for self-influence
idx = dist>0;
M = ones(size(dist));

% Matrix of cell reading
M(idx) = sinh(Rcell)./(a0*dist(idx)).*exp(Rcell-a0*dist(idx));

% Concentration in each cell
C0 = 1 + (Son-1).*cells;

% Reading of each cell
Y = M*C0;

dK = normrnd(0, noise, size(Y)); % Gaussian noise

X = (2*cells-1);
mom = sum(X.*(Y-K));

I = moranI(cells, a0*dist);

cells_out = Y > K + dK;
changed = ~isequal(cells_out, cells);



        