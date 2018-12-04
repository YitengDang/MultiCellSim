function update_cell_figure_external(h_cells, h_borders, cells, t, disp_mol, pos)  
    % plot title
    title(sprintf('t=%d', t), 'FontSize', 20);
    
    % update cell colours and positions
    if disp_mol==12
        % no signal -> white
        % signal 1 -> red
        % signal 2 -> blue 
        % signals 1&2 -> black
        
        % --update cells--
        % Update cell states
        clrs = 1-cells;
        c_all = zeros(size(cells, 1), 3); 
        c_all(:, 3) = clrs(:, 1); % signal 1 present -> Turn on blue channel
        c_all(:, 1) = clrs(:, 2); % signal 2 present -> Turn on red channel
        c_all(:, 2) = clrs(:, 1) & clrs(:, 2); % signal 2 present -> Turn on green channel
        set(h_cells, 'cdata', c_all);
        
        % Update cell positions
        set(h_cells, 'xdata', pos(:, 1), 'ydata', pos(:, 2) );
        set(h_borders, 'xdata', pos(:, 1), 'ydata', pos(:, 2) );
    else
        cells = cells(:, disp_mol);
        
        % Update cell colors
        %cells = ones(N,1);
        c_all = repmat(1-cells, 1, 3);
        set(h_cells, 'cdata', c_all);
        
        % Update cell positions
        set(h_cells, 'xdata', pos(:, 1), 'ydata', pos(:, 2) );
        set(h_borders, 'xdata', pos(:, 1), 'ydata', pos(:, 2) );
    end
end