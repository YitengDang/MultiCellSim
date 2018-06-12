function [dist, pos] = init_dist_hex(gridsizex, gridsizey)
% Try to load a previously calculated matrix. If not found calculate and
% save a new one initializing all cells in an hexagonal grid.

% Preiously saved matrices are saved in ./data/dist_matrix_hex with the
% filename <gridsizex><gridsizey>.mat
filename = fullfile(pwd, 'data', 'dist_matrix_hex', ...
    sprintf('%d%d.mat', gridsizex, gridsizey));

if exist(filename, 'file') == 2
    %disp('exist');
    tmp = load(filename);
    dist = tmp.dist;
    pos = tmp.pos;
else
    %disp('not exist');
    [pos,ex,ey] = init_cellpos_hex(gridsizex,gridsizey);
    dist = dist_mat(pos,gridsizex,gridsizey,ex,ey);
    [out, msg] = mkdir('data', 'dist_matrix_hex');
    if out~=1 % print message if output unexpected
        disp(msg);
        disp('Error msg. 1: see reference on Wiki');
    else 
        disp('Successfully created folder');
    end
    save(filename, 'pos', 'dist')
end