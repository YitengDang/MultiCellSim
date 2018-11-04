function h = plot_color_map_2D(ngenes, rows, cols)
    if nargin==1
        cols = 500;
        if ngenes==1
            rows = 75;
        else
            rows = 500;
        end
    elseif nargin==2
        cols = 500;
    end
    
    h = figure; %(11);
    h.Name = 'plot_color_map_2D';
    switch ngenes
        case 1
            % Plot 1D color map
            
           
            image=zeros(rows, cols, 3); %initialize
            for row=1:rows
                for col=1:cols
                    c = 1-(col-1)/(cols-1);
                    image(row, col, :) = [c c c];
                end
            end
            imshow(image)
            axis on
            xlabel('Gene expression');
            set(gca, 'XTick', 0:cols/5:cols, 'XTickLabel', {'0',...
                '0.2', '0.4', '0.6', '0.8', '1'} );
            set(gca, 'YTick', []);    
            set(h, 'Position', [100 100 600 200]);
        case 2 
            % Plot 2D color map
            image=zeros(rows, cols, 3); %initialize
            for row=1:rows
                for col=1:cols
                    y = 1-(row-1)/(rows-1);
                    b = 1-(col-1)/(cols-1);
                    image(row, col, :) = [y y b];
                end
            end

            imshow(image)
            xlabel('Gene 1');
            ylabel('Gene 2');
            axis on
            set(gca, 'YDir', 'normal')
            set(gca, 'XTick', 0:cols/5:cols, 'XTickLabel', {'0',...
                '0.2', '0.4', '0.6', '0.8', '1'} );
            set(gca, 'YTick', 0:rows/5:rows, 'YTickLabel', {'0',...
                '0.2', '0.4', '0.6', '0.8', '1'} );
            
            set(h, 'Position', [100 100 650 650]);
            
    end
    
    % common features
    title('Color code');
    set(gca, 'FontSize', 16);
    
end