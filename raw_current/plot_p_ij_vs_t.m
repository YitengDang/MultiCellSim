function msg = plot_p_ij_vs_t(cells_hist, fig_pos)
%{
clear variables
close all

N=64;
cells_hist = {};
for t=1:20
    cells = zeros(N, 2);
    cells(randperm(N, 34), 1) = 1;
    cells(randperm(N, 24), 2) = 1;
    cells_hist{end+1} = cells;
end
%}
%%

if isempty(cells_hist)
    msg = ' Unable to plot p(t); ';
    return
end

s = size(cells_hist{1}, 2);

if s==1 % for 1 gene: just plot p(t)
    msg = plot_p_vs_t(cells_hist, fig_pos);
    return
end

N = size(cells_hist{1}, 1);
tmax = numel(cells_hist)-1;

N_ij = zeros(tmax+1, s^2);
for t=1:tmax+1
    cells = cells_hist{t};
    base = (2.^(s-1:-1:0))';
    cells_idx = 1+cells*base;
    for state=1:s^2
        N_ij(t, state) = sum(cells_idx==state);
    end
end
%%
%
if s==1 %vectorize
    %legend_text = {'0', '1'};
    legend_text = sprintfc('%d', [0 1]);
elseif s==2
    %legend_text = {'(0,0)', '(0,1)', '(1,0)', '(1,1)'};
    legend_text = sprintfc("(%d, %d)", [0 0; 0 1; 1 0; 1 1]);
end
%}
%
h1 = figure(1);
cla(h1, 'reset');
hold on

set(gca, 'Color', [0.8 0.8 0.8]);
plot_clrs = [1 1 1; 
            0 0 1
            1 1 0
            0 0 0];
            
if length(cells_hist) < 100
    ps = 'o-'; lw = 1.5;
elseif length(cells_hist) < 500
    ps = '.-'; lw = 1;
else
    ps = '.-'; lw = 0.5;
end

for i=1:s^2
    clr = plot_clrs(i, :);
    plot(0:tmax, N_ij(:,i)/N, ps, 'LineWidth', lw, 'Color', clr);
end
%legend_text = {'a','b','c','d'};
legend(legend_text, 'Location', 'eastoutside', 'FontSize', 12);
xlabel('$$t$$', 'Interpreter', 'latex');
ylabel('$$p$$', 'Interpreter', 'latex');

set(gca, 'FontSize', 24);
xlim([0 tmax])
ylim([0 1]);

set(h1, 'Units', 'inches', 'Position', fig_pos);
%}
msg = 'Successfully plotted p_ij(t); '; 

end