function [A, phase, h] = plot_state_diagram_onecell(M_int, Con, Coff, K)
% Plots the state diagram of a single cell with ONE or TWO types of
% signalling molecules

    if all(size(M_int)==[1 1]) % One type of signalling molecule
        % Calculate A
        A = zeros(2);
        switch M_int
            case 1 % activation
                if Coff < K
                    A(1,1) = 1;
                else
                    A(1,2) = 1; % (never happens)
                end

                if Con > K
                    A(2,2) = 1;
                else
                    A(2,1) = 1;
                end
            case -1 % deactivation
                if Coff > K
                    A(1,1) = 1; % (never happens)
                else
                    A(1,2) = 1;
                end

                if Con < K
                    A(2,2) = 1;
                else
                    A(2,1) = 1;
                end
        end
        
        %% Draw state diagram (N.B. same code as in plot_state_diagram_multicell_one_signal.m)
        h = figure; %(10);
        h.Name = 'plot_state_diagram_singlecell';
        hold on
        s = [0 1];
        t = [0 0];
        Gs = digraph(A);
        nLabels = {};
        g=plot(Gs, 'XData', s, 'YData', t, 'ArrowSize', 20, 'EdgeAlpha', 1, ...
            'LineWidth', 3, 'EdgeColor', 'k', ...
            'Marker', 'o', 'MarkerSize', 100, 'NodeColor', [0.2 0.2 0.2], 'NodeLabel', nLabels);
        % Make edges dashed if state has two outgoing edges
        for i=1:2
            if sum(A(i,:))==2
                idx = find(A(i,:));
                highlight(g, i, idx, 'LineStyle', '--', 'LineWidth', 2);
            end
        end
        text(s-0.03,t,{'0', '1'}, 'Color', 'w', 'FontSize', 32)
        ax = gca;
        if A(1,2)==0 && A(2,1)==0
            axis([-0.4 1.4 -0.03 0.03]);
        else
            axis([-0.4 1.4 -0.25 0.25]);
        end
        ax.Visible = 'off';
        h.Color = [1 1 1];
        %set(ax, 'Units', 'Inches', 'Position', [0 0 9 8]);
        %set(h2, 'Units', 'Inches', 'Position', [1 1 9 8]);
        set(ax, 'Units', 'Inches', 'Position', [0 0 7 4]);
        set(h, 'Units', 'Inches', 'Position', [0.2 0.2 7 4]);
        
        %
        phase = (Con > K);
    elseif all(size(M_int)==[2 2])
        % state transitions
        X_out_sub = cell(2);
        X_out_ind = zeros(4, 1);
        A = zeros(4); % graph adjacency matrix
        for X1=0:1
            for X2=0:1 
                ind = sub2ind([2 2], X1+1, X2+1);
                %disp(ind);
                X = [X1 X2];
                C = (Con-Coff).*X+Coff;
                out = (([C; C] - K).*M_int > 0) + (1 - abs(M_int));
                X_out = prod(out, 2);
                % if no connections to a gene, output = input (remains constant)
                X_out(sum(abs(M_int), 2)==0) = X(sum(abs(M_int), 2)==0);

                X_out_2 = sub2ind([2 2], X_out(1)+1, X_out(2)+1);

                X_out_sub{X1+1,X2+1} = X_out;
                X_out_ind(ind) = X_out_2;

                % store as adjacency matrix
                A(ind, X_out_2) = 1;
            end
        end
        % Calculate phase
        phase = (repmat(Con, 2, 1) > K);

        %% Draw state diagram as directed graph
        h = figure; %(8);
        h.Name = 'plot_state_diagram_singlecell';
        hold on
        s = [0 1 0 1];
        t = [0 0 1 1];
        %A = ones(4);
        Gs = digraph(A);
        nLabels = {};
        plot(Gs, 'XData', s, 'YData', t, 'ArrowSize', 20, 'EdgeAlpha', 1, ...
            'LineWidth', 3, 'EdgeColor', 'k',...
            'Marker', 'o', 'MarkerSize', 100, 'NodeColor', [0.2 0.2 0.2], 'NodeLabel', nLabels);
        text(s-0.11,t+0.019,{'(0,0)','(1,0)','(0,1)','(1,1)'}, 'Color', 'w', 'FontSize', 32)
        ax = gca;
        axis([-0.4 1.4 -0.4 1.4]);
        ax.Visible = 'off';
        h.Color = [1 1 1];
        set(ax, 'Units', 'Inches', 'Position', [0 0 7 6]);
        set(h, 'Units', 'Inches', 'Position', [0.3 0.3 7 6]);
    end

end