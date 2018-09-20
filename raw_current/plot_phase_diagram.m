function h = plot_phase_diagram(gz, a0, rcell, M_int, K_in, Con_in, lambda12)
    set(0, 'defaulttextinterpreter', 'latex');
    
    % Phases (after readjustment)
    % 1: all ON(+) / OFF(-) 
    % 2: all OFF(+) / ON(-)
    % 3: ON->ON/activation (+) / ON-> OFF (-) 
    % 4: OFF->OFF/deactivation (+) / OFF->ON (-)
    % 5: none (activation-deactivation)
    % 6: autonomy (+) / autonomous oscillations (-)

    l = size(M_int, 1); %number of molecules
    % plots phase diagrams for each interaction of the system
    Rcell = rcell*a0;
    lambda = [1 lambda12];
    
    % Parameter range map
    num_points = 1000;
    if max(Con_in) > 40
        Con_max = ceil(max(Con_in)/10)*10;
    else
        Con_max = 40;
    end
    if max(K_in(:)) > 40
        K_max = ceil(max(K_in(:))/10)*10;
    else
        K_max = 40;
    end
    Con_vec = linspace(1, Con_max, num_points);
    K_vec = linspace(1, K_max, num_points);
    [K, Con] = meshgrid(K_vec, Con_vec);

    % calculate fN
    [dist, ~] = init_dist_hex(gz, gz);
    dist_vec = a0*dist(1,:);
    r = dist_vec(dist_vec>0); % exclude self influence
    for i=1:l % calculate signaling strength
        fN(i) = sum(sinh(Rcell)*sum(exp((Rcell-r)./lambda(i)).*(lambda(i)./r)) ); 
    end
    
    h = figure; %(9);
    h.Name = 'plot_phase_diagram';
    for idx=1:l
        subplot(l,1,idx);
        hold on

        % Make 4 limiting regions as boolean matrices
        R1 = (1+fN(idx) - K) > 0; % Everything ON
        R2 = ((1+fN(idx))*Con - K ) < 0; % Everything OFF
        R3 = ((Con + fN(idx) - K) > 0 & (1+fN(idx) - K) < 0); % ON remains ON & not all ON
        R4 = ((1+ fN(idx)*Con - K) < 0 & ((1+fN(idx))*Con - K ) > 0) ; % OFF remains OFF & not all OFF
        %R3 = (Con > K & (K-1)./Con > fN); % autonomous cells for Son > K
        %R4 = (Con <= K & K - Con < fN & (K-1)./Con > fN); % autonomous cells for Son < K

        out = R1 + 2*R2 + 3*R3 + 4*R4; % only regions 3 and 4 can overlap
        if ~isempty(find(unique(out)==0, 1))
            map_idx = 5; % activation-deactivation
            out(out==0) = 5; 
            %phase = 'non-A';
            phase = 'U';
        elseif ~isempty(find(unique(out)==7, 1)) 
            map_idx = 6; % autonomy
            out(out==7) = 5; % ON remains ON & OFF remains OFF
            phase = 'A01';
        end

        him = imagesc(K_vec, Con_vec, out);
        %set(him, 'AlphaData', out > 0); % invisible if not from any region
        % R1 -> black
        % R2 -> white
        % R3 -> green
        % R4 -> red
        % activation-deactivation -> magenta
        % autonomy -> gray
        map = [0, 0, 0
            1, 1, 1
            0, 1, 0
            1, 0, 0
            1, 1, 0
            0.5, 0.5, 0.5];
        tmp = map([1:4 map_idx], :);
        colormap(subplot(l,1,idx), tmp);
        c=colorbar;
        set(c, 'YTick', 1+2/5+4/5*(0:4));
        
        %phase_labels = {'all>K','all<K','A1','A0',phase};
        phase_labels = {'P1','P0','A1','A0',phase};
        set(c, 'TickLabels', phase_labels);
        
        %}

        % Plot current parameters as points
        markers = {'^', 'o'};
        for i1=1:l % gene influenced
            im = (M_int(i1, idx)==-1)*1 + (M_int(i1, idx)==1)*2;
            if im~=0
                scatter(K_in(i1, idx), Con_in(idx), 100, markers{im}, 'filled', 'MarkerFaceColor', [1 0 1]);
                text(K_in(i1, idx)+1, Con_in(idx)+1, num2str(i1), 'FontSize', 20, 'Color', [1 0 1]);
            end
        end
        %plot(K(idx, 2), Con(2), 'x', 'Color', [1 0 1]);
        %}
        % Plot lines
        %{
        hold on
        Con2 = linspace(1,30,1000);
        K2 = linspace(1,20,1000);
        plot([fN+1 fN+1], [1 Con2(end)], 'k--', 'LineWidth', 1.5) % all ON region
        plot([fN+1 fN+1], [1 Con2(end)], 'g--', 'LineWidth', 1.5) % all ON region
        plot(K2, K2./(1+fN), 'b--', 'LineWidth', 1.5) % all OFF region
        plot(K2, K2-fN, 'Color', [247 145 52]/256, 'LineStyle', '--', 'LineWidth', 1.5) % autonomous 1
        plot(K2, (K2-1)/fN, 'r--', 'LineWidth', 1.5) % autonomous 2
        % plot(K, 2*K/(1+fN)-1, '--k', 'LineWidth', 1.5) % (Con+1)(1+fN) = 2K
        hold off
        %}

        % adjust the graph
        set(gca,'ydir', 'normal', 'FontSize', 24)
        xlabel('K', 'FontSize', 24)
        ylabel('$$C_{ON}$$', 'FontSize', 24)
        ylim([1 Con_max])
        xlim([1 K_max])
        %title(sprintf('$$f_N = %.3f, a_0 = %.2f$$', fN, a0), 'FontSize', 30)
        %title(sprintf('Interaction $$%d \\leftarrow %d$$', idx, idx2))
        title(sprintf('Molecule %d', idx));

        % set ticks
        if K_max<100 && Con_max<100
            tick_K = [1 10:10:K_max];
            tick_Con = [1 10:10:Con_max];
        else
            tick_K = [1 K_max/10:K_max/10:K_max];
            tick_Con = [1 Con_max/10:Con_max/10:Con_max];
        end
        set(gca, 'xtick', tick_K, 'ytick', tick_Con)

        
    end

    set(h, 'Units', 'Inches', 'Position', [1 1 10 9]);
    
end
%{ 
% Code for testing
clear all
close all

% input parameters
gz = 15;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;
M_int = [1 -1; -1 1];
K_in = [12 4; 20 5];
Con_in = [24 20];
lambda = [1 2];

plot_phase_diagram_1(gz, a0, rcell, M_int, K_in, Con_in, lambda(2));
%}
