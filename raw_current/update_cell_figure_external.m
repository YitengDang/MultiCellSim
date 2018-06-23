function update_cell_figure_external(hin, pos, cells_in, cell_type, t, disp_mol)
    % Plot the cell state in a diagram according to their type. If
    % cell_type(i) = 0, the cell is plotted as circle, if 1 it is plotted
    % as square. hin is the figure handle.
    
    % clear figure and hold on to keep the graphics
    %clf(hin,'reset');
    figure(hin);
    title(sprintf('Time = %d', t), 'FontSize', 24);
    set(gca,'YTick',[],'XTick',[]);
    box on
    hold on
    %title(strcat('$$\langle X_i \rangle$$',...
    %    sprintf('= %.2f, Time: %d', mean(cells), t)), ...
    %    'FontSize', 24, 'Interpreter', 'latex');
    %% calculate figure dimensions
    N = size(pos,1);
    gz = sqrt(N);
    Lx = 560;
    Ly = (sqrt(3)/2*(gz-1)+2)/(3/2*(gz+1))*Lx;
    a0 = Lx/(3/2*(gz+1));
    Rcell = 0.5*a0;
    %%
    % set image properties
    %h = gcf;
    %hold on
    xlim([-1 pos(end,1)+1]);
    ylim([-1 pos(end,2)+1]);
    set(hin, 'Units', 'points');
    set(hin, 'Position', [50 50 Lx Ly]);
    set(gca, 'color', [0.8 0.8 0.8]) % background color
     
    % get the right cells to display
    if disp_mol==12
        cells = cells_in;
        clrs = 1-cells;
        c_all = zeros(N, 3); 
        c_all(:, 3) = clrs(:, 1); % signal 1 present -> Turn on blue channel
        c_all(:, 2) = clrs(:, 2); % signal 2 present -> Turn on green channel
        c_all(:, 1) = clrs(:, 2); % signal 2 present -> Turn on red channel
    else
        cells = cells_in(:, disp_mol);
        c_all = repmat(1-cells, 1, 3); % all colours
        %markers = {'o', 's'};
    end
    
    % plot cells
    % (1) plot with imagesc
    %cells_2 = reshape(1-cells, gz, gz);
    %imagesc(pos(:, 1), pos(:, 2), cells_2);
    
    % (2) plot as loop
    %for i = 1:N
    %    plot(pos(i,1), pos(i,2), 'Marker', markers{cell_type(i)+1}, 'MarkerFaceColor',...
    %        c_all(i, :), 'MarkerEdgeColor', 'k', 'MarkerSize', Rcell); 
    %end
    
    % (3) use scatter
    %cell_type = zeros(N, 1);
    %cell_type(randperm(N, round(0.5*N))) = 1;
    
    clr_k = zeros(N, 3);
    idx0 = (cell_type==0);
    p0 = scatter(pos(idx0,1), pos(idx0,2), Rcell^2, c_all(idx0, :), 'filled', 'o');
    scatter(pos(idx0,1), pos(idx0,2), Rcell^2, clr_k(idx0, :), 'o'); % plot cell boundaries
    %{
    idx1 = (cell_type==1);
    p1 = scatter(pos(idx1,1), pos(idx1,2), Rcell^2, c_all(idx1, :), 'filled', 's');
    scatter(pos(idx1,1), pos(idx1,2), Rcell^2, clr_k(idx1, :), 's'); % plot cell boundaries
    %}
    hold off
    % Add color bar for Xi
    c = colorbar();
    c.Ticks = 0:0.2:1;
    set(c, 'FontSize', 14);
    ylabel(c, '$$X_i$$', 'Interpreter', 'latex', 'FontSize', 16,...
        'Rotation', 90);
    map = repmat((1:-0.01:0)', 1, 3);
    %map = 'gray';
    colormap(map);
    %}