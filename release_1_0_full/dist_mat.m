function dd = dist_mat(pos, Lx, Ly, ex, ey)
% Calculate the distance matrix. Lx and Ly are set in case of periodic
% boundary (box size). If there is a connected pattern, then you should also
% provide the unit vectors and Lx and Ly are interpreted as the grid sizes
% For hexagonal lattice ex = [1 0] and ey = 0.5*[1 sqrt(3)]

% pos is the positioning of all the cells, attributed before this function
% This function simply calculate a matrix with the distance between cells.

% No periodic boundary in case L's are not given.
if nargin == 1
    Lx = 0;
    Ly = 0;
    ex = [];
    ey = [];
else
    if nargin == 3
        ex = [];
        ey = [];
    else
        if nargin ~= 5
            sprintf('Wrong number of parameters')
        end
    end
end

% get the number of cells
n = size(pos,1);
dd = zeros(n);

for i = 1 : n
    for j = i+1 : n
        diff_pos = pos(i,:) - pos(j,:);
        % Apply periodic boundary. If difference in the given direction is
        % bigger than half of the given dimension, then the round function
        % leads to 1 and we suctract the size of the box to generate the
        % periodicity. The same argument holds for negative differences.
        if Lx ~= 0 && Ly ~= 0
            if isempty(ex) && Lx ~= 0
                diff_pos(1) = diff_pos(1) - round(diff_pos(1)./Lx).*Lx;
            end
            if isempty(ey) && Ly ~= 0
                diff_pos(2) = diff_pos(2) - round(diff_pos(2)./Ly).*Ly;
            end
            % If we have a generated pattern, then we need to project the
            % distances into the unit vectors and compare to the grid sizes
            % 
            if ~(isempty(ex) && isempty(ex))
                % Find the distances in terms of the unit vectors. In the
                % end we have diff_pos = projx*ex+projy*ey.
                exey = dot(ex,ey);
                if exey==0
                    projx = dot(diff_pos,ex);
                    projy = dot(diff_pos,ey);
                else
                    eperp = ey - ex/exey;
                    projx = dot(diff_pos,eperp)/dot(ex,eperp);
                    projy = dot(diff_pos,ey) - projx*exey;
                end
                % Now we can correct the lengths by comparing the
                % projections with the gridsizes
                diff_pos = diff_pos - round(projx/Lx).*Lx*ex;
                diff_pos = diff_pos - round(projy/Ly).*Ly*ey;
            end
        end

        dd(i,j) = sqrt(sum(diff_pos.^2));
        dd(j,i) = dd(i,j);
    end
end
