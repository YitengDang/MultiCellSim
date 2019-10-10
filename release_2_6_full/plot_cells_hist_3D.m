% Plots cells_hist as a 3D plot
N = size(pos, 1);

figure;
hold on
for i=1:1
    cells = cells_hist{i};
    
    % get cell colours
    clrs = 1-cells;
    c_all = zeros(size(cells, 1), 3); 
    c_all(:, 3) = clrs(:, 1); % signal 1 present -> Turn on blue channel
    c_all(:, 2) = clrs(:, 2); % signal 2 present -> Turn on green channel
    c_all(:, 1) = clrs(:, 2); % signal 2 present -> Turn on red channel
    
    % plot all
    %scatter3(pos(:,1), pos(:,2), i*ones(N, 1), 50, c_all, 'filled', 'MarkerFaceAlpha', 0.5);
    
    % plot only one type of cell
    idx = find( all(cells==[0 0], 2) );
    scatter3(pos(idx,1), pos(idx,2), i*ones(numel(idx), 1), 150, c_all(idx,:), 'filled', 'MarkerFaceAlpha', 0.5);

end
view([45 30]);
set(gca, 'Color', [0.8 0.8 0.8]);

%%
