function [I, theta] = moranI(cells, dist)
% This function calculates the modified Moran I or spatial order parameter.
% The weights are assumed to be the interaction f_ij

cells_pm = 2*cells - 1; % ON cells: 1, OFF cells: -1
cell_mean = mean(cells_pm); % <Xi>

% meshgrid replicates the cells states in x and y directions
[cells_matx, cells_maty] = meshgrid(cells_pm, cells_pm);

% For the Moran I we use the factor of the diffusion equation as weight
% Account for self-influence
idx = dist>0;
M = zeros(size(dist));

% Matrix of cell reading
M(idx) = exp(-dist(idx))./dist(idx);

% Sum of all weights
w_sum = sum(sum(M));
theta = sum(sum(M.*cells_matx.*cells_maty))/w_sum;

% Multiply every weight by the cell state in a crossed fashion
tmp = sum(sum(M.*(cells_matx-cell_mean).*(cells_maty-cell_mean)));
if tmp==0
    I = 0;
else
    % cells state variance (using N as normalization)
    cells_var = var(cells_pm,1);
    I = sum(sum(tmp))/w_sum/cells_var;
end