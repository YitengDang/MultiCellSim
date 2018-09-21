function [heat_matrix, msg, h] = plot_time_state_heatmap(cells_hist, genes)
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
Gridsize = sqrt(N);
heat_matrix = zeros(length(cells_hist),N);

if genes == 1 || genes == 2
    for t = 1:length(cells_hist)
        h_index = 1;
        for i = Gridsize:-1:1
            for j = i:Gridsize:Gridsize*Gridsize-(Gridsize - i)
                 heat_matrix(t,h_index) = cells_hist{t}(j,genes);
                 h_index = h_index + 1;
            end
        end
    end
elseif genes == 12
    for t = 1:length(cells_hist)
        h_index = 1;
        for i = Gridsize:-1:1
            for j = i:Gridsize:Gridsize*Gridsize-(Gridsize - i)
                cell_state = cells_hist{t}(j,:);
                %disp(cell_state);
                heat_matrix(t,h_index) = sum(cell_state.*[1 2]);
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
% White  -> 00
% yellow -> 10
% blue   -> 01
% black  -> 11
h = figure;
h.Name = 'plot_kymograph';
imagesc(heat_matrix)
xlabel('Cell')
ylabel('Time')
if genes==1 || genes==2
    cmap = flipud(gray);
    colormap(cmap)
    title(sprintf('Gene %d',genes))
elseif genes==12
    mymap = [1 1 1
    1 1 0
    0 0 1
    0 0 0];
    colormap(mymap)
    title('Both genes');
end
%colorbar
set(gca, 'FontSize', 16);

msg = 'Successfully plotted kymograph';
end