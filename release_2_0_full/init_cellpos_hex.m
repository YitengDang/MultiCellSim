% Initialize the cell positions in hexagonal lattice with specified grid
% pos is a (N,2) array with the (x,y) coordinates of all cells.

function [pos,e1,e2] = init_cellpos_hex(gridsize_x, gridsize_y)
    
    e1 = [1/2 sqrt(3)/2];
    e2 = [1 0];
    
    pos = zeros(gridsize_x*gridsize_y, 2);
    for i = 1:gridsize_x
        for j = 1:gridsize_y
            pos((i-1)*gridsize_y+j,:) = e1*(j-1) + e2*(i-1);
        end
    end
    