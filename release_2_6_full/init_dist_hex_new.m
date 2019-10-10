function [dist, pos] = init_dist_hex_new(gridsize, mcsteps, rcell, periodic_bc)
% Try to load a previously calculated matrix. If not found calculate and
% save a new one initializing all cells in an hexagonal grid.
% 
% Preiously saved matrices are saved in ./data/dist_matrix_hex with the
% filename <gridsizex><gridsizey>.mat
if all(periodic_bc)
    fname_str = sprintf('%d%d.mat', gridsize, gridsize);
else
    fname_str = sprintf('%d%d_periodic_%d_%d.mat', gridsize, gridsize,...
        periodic_bc(1), periodic_bc(2));
end
filename = fullfile(pwd, 'data', 'dist_matrix_hex', fname_str);

if exist(filename, 'file') == 2
    %disp('exist');
    tmp = load(filename);
    dist = tmp.dist;
    pos = tmp.pos;
else
    %disp('not exist');
    [pos, dist, ~, ~] = ...
        initial_cells_random_markov(gridsize, mcsteps, rcell, periodic_bc);
    [out, msg] = mkdir('data', 'dist_matrix_hex');
    if out~=1 % print message if output unexpected
        disp(msg);
        disp('Error: could not make folder ./data/dist_matrix_hex');
    else 
        disp('Successfully created folder');
    end
    save(filename, 'pos', 'dist')
end