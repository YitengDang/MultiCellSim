function [heat_matrix, msg, h] = plot_kymograph(cells_hist, genes)
% This function creates a heatmap of the cellstates. On the y-axis are all
% the simulated time points, on the x-axis are all the cells. It should be
% noted that the cells of the heatmap from left to right,
% represent the grid cells from upper left of the grid (as displayed in the app)
% to lower right of the grid. In the heatmap all gridrows are placed
% next to each other. In oder to acchieve this, this function cannot just plot the
% cell states in order as they appear in cell_hist. In cell_hist the cells
% are numbered from the lowerleft upwards.
% This function can plot the state of one gene (black = on, white = off) or
% of both genes together (colourscheme below).
% Cells_hist is an array in which each element consists of a vector holding
% all cell states. 
% If gene is 1 or 2, the heatmap shows only gene 1 or 2 in black and white. If gene is 12 a
% colour scheme is used. In case of a single gene 1 should be entered.
msg = '';
h = [];

N = size(cells_hist{1},1);
gz = sqrt(N);

if genes == 1 || genes == 2
    heat_matrix = zeros(length(cells_hist), N);

    for t = 1:length(cells_hist)
        h_index = 1;
        for i = gz:-1:1
            for j = i:gz:gz*gz-(gz - i)
                 heat_matrix(t,h_index) = cells_hist{t}(j,genes);
                 h_index = h_index + 1;
            end
        end
    end
elseif genes == 12
    heat_matrix = zeros(length(cells_hist), N, 3); % RGB colors

    for t = 1:length(cells_hist)
        h_index = 1;
        
        % get colors for cells
        cells = cells_hist{t};
        clrs = 1-cells;
        c_all = zeros(size(cells, 1), 3); 
        c_all(:, 3) = clrs(:, 1); % signal 1 present -> Turn on blue channel
        c_all(:, 1) = clrs(:, 2); % signal 2 present -> Turn on red channel
        
        c_all(:, 2) = clrs(:, 1) & clrs(:, 2); % both signals present -> Turn on green channel
        
        for i = gz:-1:1
            for j = i:gz:gz*gz-(gz - i)
                %cell_state = cells(j,:);
                %disp(cell_state);
                % heat_matrix(t, h_index) = sum(cell_state.*[1 2]); % imagesc version
                
                heat_matrix(t, h_index, :) = c_all(j, :); % imagesc version
                
                %{
                gene_1 = cells_hist{t}(j,1); 
                gene_2 = cells_hist{t}(j,2);
                sum_genes = gene_1 + gene_2;
            
                if sum_genes == 0 % both off
                    heat_matrix(t,h_index) = 0;
                elseif sum_genes == 2 % both on
                    heat_matrix(t,h_index) = 3;
                elseif gene_1 == 1 % only gene 1 on
                    heat_matrix(t,h_index) = 1;
                elseif gene_2 == 1 % only gene 2 on
                    heat_matrix(t,h_index) = 2;
                else
                    error('Error in determining state of the two signalling genes.')
                end
                %}
                h_index = h_index + 1;
            end
        end
    end
else
    error('The genes should be entered as either: 1 (only gene 1 is plotted), 2 (only gene 2 is plotted), or 12 (both plotted)')
end
%% Plot 
if genes==1 || genes==2
    h = figure;
    h.Name = 'plot_kymograph';
    imagesc(heat_matrix)
    xlabel('Cell number')
    ylabel('Time')

    cmap = flipud(gray);
    colormap(cmap)
    title(sprintf('Gene %d',genes))
    
elseif genes==12
    % White  -> 00
    % red    -> 10
    % blue   -> 01
    % black  -> 11
    h = figure;
    h.Name = 'plot_kymograph';
    image(heat_matrix)
    xlabel('Cell number')
    ylabel('Time')
    
    hold on
    for ii=1:gz-1
        plot(ii*[gz gz]+0.5, [0.5 size(heat_matrix, 1)+0.5], '--', 'Color', [0.5 0.5 0.5]);
    end
    
    cmap = flipud(gray);
    colormap(cmap)
    title(sprintf('Gene %d',genes))
    
    if length(cells_hist)<10
       ystep = 1;
    elseif length(cells_hist)<50
       ystep = 5;
    elseif length(cells_hist)<100
       ystep = 10;
    elseif length(cells_hist)<500
       ystep = 50;
    elseif length(cells_hist)<1000
       ystep = 100;
    else
       ystep = 500;
    end
    
    set(gca, 'XTick', 2*gz:2*gz:N);
    set(gca, 'YTick', 1:ystep:size(heat_matrix, 1),...
        'YTickLabel', sprintfc('%d', 0:ystep:size(heat_matrix, 1)-1) );
    %{
    mymap = [1 1 1
    1 1 0
    0 0 1
    0 0 0];
    colormap(mymap)
    %}
    title('Both genes');
end
%colorbar
set(gca, 'FontSize', 16);
%set(h, 'Units', 'Inches', 'Position', [0 0 10 8]);

msg = 'Successfully plotted kymograph';
end