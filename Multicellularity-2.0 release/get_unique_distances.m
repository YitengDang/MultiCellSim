function dist_val = get_unique_distances(dist, eps)
% Get all the possible distances in the distance matrix of the distribution
% of cells in space

% eps is the acceptable difference when determining if two distances are
% equal or not
    if nargin < 2
        eps = 1e-4;
    end

    dists = unique(dist);
    dist_val = [];
    for i = 1:numel(dists)
        if isempty(find(dist_val < dists(i)+eps & dist_val > dists(i)-eps,1))
            dist_val(end+1) = dists(i);
        end
    end
end