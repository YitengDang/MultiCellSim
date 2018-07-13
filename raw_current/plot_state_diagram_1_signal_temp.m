%Phases:
% 1: none (activation-deactivation)
% 2: all ON(+) / OFF(-)
% 3: ON->ON/activation (+) / ON-> OFF (-) 
% 4: OFF->OFF/deactivation (+) / OFF->ON (-)
% 5: all OFF(+) / ON(-)
% 6: autonomy (+) / autonomous oscillations (-)
phase = 2;
M_int = -1;
idx = (M_int==1) + (M_int==-1)*2;

A_map = cell(2, 6);
% activation 
A_map{1,1} = ones(2);
A_map{1,2} = [0 1; 0 1];
A_map{1,3} = [1 1; 0 1];
A_map{1,4} = [1 0; 1 1];
A_map{1,5} = [1 0; 1 0];
A_map{1,6} = [1 0; 0 1];
% repression
A_map{2,1} = ones(2);
A_map{2,2} = [1 0; 1 0];
A_map{2,3} = [1 1; 1 0];
A_map{2,4} = [0 1; 1 1];
A_map{2,5} = [0 1; 0 1];
A_map{2,6} = [0 1; 1 0];

A = A_map{idx, phase};
%% Draw state diagram
h10 = figure(10);
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
h10.Color = [1 1 1];
%set(ax, 'Units', 'Inches', 'Position', [0 0 9 8]);
%set(h2, 'Units', 'Inches', 'Position', [1 1 9 8]);
set(ax, 'Units', 'Inches', 'Position', [0 0 7 4]);
set(h10, 'Units', 'Inches', 'Position', [0.2 0.2 7 4]);
