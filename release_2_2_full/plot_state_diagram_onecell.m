function plot_state_diagram_onecell(M_int, Con, Coff, K)
    
if numel(M_int)==1
    return
end
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

%% Draw state diagram as directed graph
h8 = figure(8);
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
h8.Color = [1 1 1];
set(ax, 'Units', 'Inches', 'Position', [0 0 7 6]);
set(h8, 'Units', 'Inches', 'Position', [0.3 0.3 7 6]);


end