function h = plot_circuit(gz, a0, rcell, M_int, Con, Coff, K, lambda12)
% plots genetic circuit
% h: figure handle

Rcell = rcell*a0;
[dist, ~] = init_dist_hex(gz, gz);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence

h = figure; %(7);
h.Name = 'plot_circuit';
hold on 
if all(size(M_int)==[1 1])
    fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength
    interactions = zeros(1,2);
    interactions(1) = (Con*(1+fN) > K);
    interactions(2) = (Coff*(1+fN) < K);
    
    G = digraph(1);
    
    color = (M_int==1)*[0 0 1] + (M_int==-1)*[1 0 0];
    g = plot(G, 'XData', 0, 'YData', 0, 'ArrowSize', 20, 'EdgeAlpha', 1, ...
        'EdgeColor', color, 'LineWidth', 3, ...
        'Marker', 'o', 'MarkerSize', 50, 'NodeColor', 'k', 'NodeLabel', {});

    if ~interactions(1,1) % always above threshold
        highlight(g, 1, 1, 'LineStyle', ':');
    elseif ~interactions(1,2) % always below threshold
        highlight(g, 1, 1, 'LineStyle', '--');
    end
    
    p1 = plot(0, 0, 'b', 'LineWidth', 2);
    p2 = plot(0, 0, 'r', 'LineWidth', 2);
    %legend([p1 p2], {'activation', 'repression'}, 'FontSize', 16, 'Location', 'best');
    
    ax = gca;
    ax.Visible = 'off';
    h.Color = [1 1 1];

    set(ax, 'Units', 'Inches', 'Position', [0 0 4 4]);
    set(h, 'Units', 'Inches', 'Position', [1 1 4 4]);

elseif all(size(M_int)==[2 2])
    
    lambda = [1 lambda12];

    % calculate fN
    fN = zeros(2,1);
    fN(1) = sum(sinh(Rcell)*sum(exp((Rcell-r)./lambda(1)).*(lambda(1)./r)) ); % calculate signaling strength
    fN(2) = sum(sinh(Rcell)*sum(exp((Rcell-r)./lambda(2)).*(lambda(2)./r)) ); % calculate signaling strength

    % conditions
    interactions = zeros(2,2,2); %0: ON=OFF, 1: ON~=OFF
    for i=1:2
        for j=1:2
            interactions(i,j,1) = (Con(j)*(1+fN(j)) > K(i,j));
            interactions(i,j,2) = (Coff(j)*(1+fN(j)) < K(i,j));
        end
    end

    %% Draw gene circuit
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
            if M_int(i,j)~=0
                if ~interactions(i,j,1) % always above threshold
                    highlight(g, j, i, 'LineStyle', ':');
                elseif ~interactions(i,j,2) % always below threshold
                    highlight(g, j, i, 'LineStyle', '--');
                end
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
    %legend([p1 p2], {'activation', 'repression'}, 'FontSize', 16, 'Location', 'best');
    
    colormap(map);
    ax.Visible = 'off';
    h.Color = [1 1 1];

    set(ax, 'Units', 'Inches', 'Position', [0 0 8 4]);
    set(h, 'Units', 'Inches', 'Position', [1 1 8 4]);

    %}
end

end