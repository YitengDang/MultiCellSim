function [cells_out, test, I_out, t_out] = ...
    generate_I_new(cells_in, I_min, I_max, dist, a0, maxsteps)
% This function tries to generate a pattern with the same fraction of ON cells as
% the vector cells_in, with spatial order I_min < I < I_max. This is
% done by trying to turn ON and OFF random cells, taking into account their neighbours
% For instance, to increase I, it tries to turn on OFF-cells that are close
% to ON cells. 
% The output is not guaranteed to produce a configuration with I in the
% correct range, but will try to get as close to the set range as possible.

% Input:
% cells - a configuration of cells with the desired fraction of ON cells
% (p)
% I_min - lower bound for allowed I
% I_max - upper bound for allowed I
% dist - matrix of distances between cells
% a0 - dimensionless distance between neighbouring cells
% maxsteps - maximum number of times to try switching cell states (to prevent
% an infinite loop)

% Output:
% cells_out - generated configuration of cells
% I_out - final value of I.
% t_out - number of trials to generate final pattern

% Suggestion: For I_max = I_min + epsilon, then epsilon = 0.01 should give
% fairly good results

%{
clear variables
close all

gridsize = 11;
N = gridsize^2;
cells_in = zeros(N, 1);
p0 = 0.2;
iniON = round(p0*N);
cells_in(randperm(N, iniON)) = 1;
I_min = 0.26;
I_max = 0.27;

% Initialize parameters
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);
%}
%%
if nargin<6
    maxsteps = 10000; 
end

% Get the number of cells
N = numel(cells_in);
[I, ~] = moranI(cells_in, a0*dist);

% Check that the lattice is not all ON/OFF or single ON/OFF
if sum(cells_in)==N || sum(cells_in)==0 || sum(cells_in)==1 || sum(cells_in)==N-1
    check = false;
else 
    check = true;
end

% Determine whether to increase or decrease I
increase = (I < I_min); % 0: decrease I, 1: increase I

% Get first neighbors
eps = 1e-5;
dist_vec = get_unique_distances(dist, eps);
dist1 = dist_vec(2);
first_nei = 1*(dist < dist1+eps & dist > dist1-eps);

%cell_type = zeros(N,1);
%h1=figure(1);
%update_cell_figure(h1, pos, 1, cells_in, cell_type, 0);
%fprintf('I = %.3f \n', I_1)

%% Increase I
%--- start while loop------
t = 0;
while (I < I_min || I > I_max) && t < maxsteps && check
    %k = waitforbuttonpress;
    t = t+1;
    cells_new = cells_in; 

    % number of neighbors that are ON
    nei_ON = first_nei*cells_new;

    % get ON cell with min./max. # ON neighbours
    nei_ON_1 = nei_ON(logical(cells_new));
    cond1 = 1;
    if increase % distinguish between increasing and decreasing I
        idx_temp = find(nei_ON_1 < 3);
        if isempty(idx_temp) % if no good cells
            idx_temp = find(nei_ON_1 == min(nei_ON_1)); % select most suitable
            cond1 = 0;
        end
    else
        idx_temp = find(nei_ON_1 > 3);
        if isempty(idx_temp) 
            idx_temp = find(nei_ON_1 == max(nei_ON_1), 1);
            cond1 = 0;
        end
    end
    idx_ON = datasample(idx_temp, 1);
    idx_1 = find(cells_new, idx_ON); %original cell index

    % get OFF cell with max./min. # ON neighbours
    nei_ON_0 = nei_ON(~logical(cells_new));
    cond2 = 1;
    if increase
        idx_temp = find(nei_ON_0 > 3);
        if isempty(idx_temp) 
            idx_temp = find(nei_ON_0 == max(nei_ON_0));
            cond2 = 0;
        end
    else
        idx_temp = find(nei_ON_0 < 3);
        if isempty(idx_temp) 
            idx_temp = find(nei_ON_0 == min(nei_ON_0));
            cond2 = 0;
        end
    end
    idx_OFF = datasample(idx_temp, 1);
    idx_0 = find(~cells_new, idx_OFF); %original cell index
    
    % ON -> OFF
    cells_new(idx_1(end)) = ~cells_new(idx_1(end));
    %h2=figure(2);
    %update_cell_figure(h2, pos, 1, cells_out, cell_type, 0);
    %fprintf('I = %.3f \n', moranI(cells_out, a0*dist))

    % OFF -> ON
    cells_new(idx_0(end)) = ~cells_new(idx_0(end));
    %h3=figure(3);
    %update_cell_figure(h3, pos, 1, cells_out, cell_type, 0);
    %fprintf('I = %.3f \n', moranI(cells_out, a0*dist))

    % if not sure, check whether I has increased/decreased
    I_new = moranI(cells_new, a0*dist);
    if cond1 && cond2 % conditions met
        cells_in = cells_new; % accept change
        I = I_new;
    elseif increase && (I_new >= I) % conds not met, but new I increased as required
        cells_in = cells_new; % accept change        
        I = I_new;
    elseif (I_new <= I) % conds not met, but new I decreased as required
        cells_in = cells_new;
        I = I_new;
    else
        %disp('rejected!');
    end
    % (else: reject change)
    increase = (I < I_min); % reanalyze whether to increase or decrease I
end

%% Output
cells_out = cells_in;
test = (I > I_min) && (I < I_max); %Check whether final I is in desired range
I_out = I;
t_out = t;
%%
%h2=figure();
%update_cell_figure(h2, pos, 1, cells_out, cell_type, 0);