function [cells_out, changed] = ...
    update_cells_noise_hill(cells, dist, Son, K, a0, Rcell, noise, hill, prec)
% Update cells using noise in a positive feedback loop with finite hill
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

%dY = noise*sqrt(Y); % noise a la Berg-Purcell
%Y = Y + dY.*(2*rand(size(Y))-1);
dY = normrnd(0, noise, size(Y)); % Gaussian noise
Y = Y+dY;

cells_out = Y.^hill./(Y.^hill + K.^hill);
changed = ~isequal(round(cells_out, prec), round(cells, prec));



        