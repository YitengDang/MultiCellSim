function plot_circuit(M_int, Con, K)
% plots genetic circuit

if all(size(M_int)==[1 1])
    % to be completed
elseif all(size(M_int)==[2 2])
    % conditions
    interactions = zeros(2); %0: ON=OFF, 1: ON~=OFF
    for i=1:2
        for j=1:2
            interactions(i,j) = (Con(j) > K(i,j));
        end
    end

    %% Draw gene circuit
    h7 = figure(7);
    hold on
    M2 = M_int'; 
    G = digraph(abs(M2));
    s = [-1 1];
    t = [0 0];
    
    edgecolors = M_int(M_int~=0);

    % (1) automatic layout
    %plot(G, 'Layout', 'circle', 'ArrowSize', 20, 'EdgeAlpha', 1, ...
    %    'EdgeCData', edgecolors, 'LineWidth', 3, 'NodeLabel', {},...
    %    'Marker', 'o', 'MarkerSize', 100, 'NodeColor', 'k');

    % (2) manual layout
    g = plot(G, 'XData', s, 'YData', t, 'ArrowSize', 20, 'EdgeAlpha', 1, ...
        'EdgeCData', edgecolors, 'LineWidth', 3, ...
        'Marker', 'o', 'MarkerSize', 100, 'NodeColor', 'k', 'NodeLabel', {});

    % Node Labels
    text(s-0.03,t+0.005,{'1', '2'}, 'Color', 'w', 'FontSize', 40)

    % Edge styles
    % solid lines -> interaction threshold reached
    % dashed lines -> interaction threshold not reached
    for i=1:2
        for j=1:2
            if ~interactions(i,j) && M_int(i,j)~=0
                highlight(g, j, i, 'LineStyle', ':')
            end
        end
    end

    % Edge Labels (abandoned)
    %{
    eLabels = {'x', 'o'};
    for i=1:2
        for j=1:2
            if M2(i,j)==0
                continue
            end
            if interactions(i,j)==0
                eLabels{end+1} = eLabels{1};
            elseif interactions(i,j)==1
                eLabels{end+1} = eLabels{2};
            end
        end
    end
    s2 = [-2.2 -0.05; -0.05 2.1];
    t2 = [0.01 -0.15; 0.15 0.01];
    text(s2(M2~=0), t2(M2~=0), eLabels, 'Color', [0 0.7 0], 'FontSize', 24); 
    %}

    ax = gca;
    %axis([-1 3 -0.2 0.2]);
    %colormap(hsv(2));
    % manually set colors for when there is only activation/repression
    if sum(unique(M_int)==1)==1&&sum(unique(M_int)==1)==1
        map = [1, 0, 0 
            0, 0, 1];
    elseif sum(unique(M_int)==1)==1
        map = [0, 0, 1];
    elseif sum(unique(M_int)==-1)==1
        map = [1, 0, 0];
    else
        map = [0, 0, 0];
    end
    
    % plot legend
    p1 = plot(0, 0, 'b', 'LineWidth', 2);
    p2 = plot(0, 0, 'r', 'LineWidth', 2);
    legend([p1 p2], {'activation', 'repression'}, 'FontSize', 16, 'Location', 'best');
    
    colormap(map);
    ax.Visible = 'off';
    h7.Color = [1 1 1];

    set(ax, 'Units', 'Inches', 'Position', [0 0 8 4]);
    set(h7, 'Units', 'Inches', 'Position', [1 1 8 4]);

    %}
end

end