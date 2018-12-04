function [h_cells, h_borders] = reset_cell_figure(ax, pos, rcell)
    % h_cells: handle of the scatter plot of cell states
    % h_borders: handle of the scatter plot of cell borders
    
    % format plot
    %ax = hin;
    cla(ax);
    
    % calculate figure dimensions
    N = size(pos,1);
    gz = sqrt(N);
    
    % all sizes in units of pixels
    Sx = 500; %ax.Position(3); %512;
    Sy = sqrt(3)/2*Sx; %(sqrt(3)/2*(gz-1)+2)/(3/2*(gz+1))*Sx;
    set(ax, 'Units', 'points', 'Position', [100 100 1.2*Sx 1.2*Sy]);
    
    % set image properties (in terms of plot units, e.g. 0 <= x <= 1) 
    set(gca, 'Units', 'points', 'Position', [0.1*Sx 0.1*Sy Sx Sy]);
    set(gca, 'YTick', [], 'XTick', [], 'Color', [0.8 0.8 0.8]);
    title(gca, 'Simulate dynamics', 'FontSize', 20);
    Lx = 1;
    Ly = sqrt(3)/2*Lx;
    %d = 1/a0*Lx/(gz); %rcell*Lx/(gz+1);
    a0_px = 1/(gz); 
    d = a0_px/2;
    xlim([-d Lx+d]);
    ylim([-d Ly+d]);
    
    % set a0, Rcell in terms of figure sizes
    h_gca = gca;
    set(h_gca, 'Units', 'points');
    Sx_px = h_gca.Position(3);
    Lx_px = (Lx/(Lx+2*d))*Sx_px;
    a0_px = Lx_px/(gz); %Lx/(3/2*(gz+1));
    Rcell_px = rcell*a0_px;
    
    %% --plot cells--
    hold on
    % colours
    c_all = ones(N, 3);
    clr_k = zeros(N, 3); % black boundaries
    %markers = {'o', 's'};

    h_cells = scatter(pos(:,1), pos(:,2), (2*Rcell_px)^2, c_all, 'filled', 'o');
    h_borders = scatter(pos(:,1), pos(:,2), (2*Rcell_px)^2, clr_k, 'o'); % plot cell boundaries
    
    % Plot box outline
    plot([0 Lx], [0 0], 'k--');
    plot([0 Lx], [Ly Ly], 'k--');
    plot([0 0], [0 Ly], 'k--');
    plot([Lx Lx], [0 Ly], 'k--');
    
    % plot round and square cells separately
    %{
    idx0 = (cell_type==0);
    p0 = scatter(ax,pos(idx0,1), pos(idx0,2), Rcell^2, c_all(idx0, :), 'filled', 'o');
    scatter(ax,pos(idx0,1), pos(idx0,2), Rcell^2, clr_k(idx0, :), 'o'); % plot cell boundaries

    idx1 = (cell_type==1);
    p1 = scatter(ax,pos(idx1,1), pos(idx1,2), Rcell^2, c_all(idx1, :), 'filled', 's');
    scatter(ax,pos(idx1,1), pos(idx1,2), Rcell^2, clr_k(idx1, :), 's'); % plot cell boundaries
    %}
    hold off
    %% Add color bar for Xi
    %{
    c = colorbar;
    c.Ticks = 0:0.2:1;
    set(c, 'FontSize', 12);
    %c.Label.String = '$$X_i$$';
    %ylabel(c, '$$X_i$$', 'Interpreter', 'latex', 'FontSize', 16,...
    %    'Rotation', 90);
    caxis([0 1]);
    map = repmat((1:-0.01:0)', 1, 3);
    %map = 'gray';
    colormap(map);
    %}
end