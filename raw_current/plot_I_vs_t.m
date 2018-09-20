function [msg, h] = plot_I_vs_t(cells_hist, a0, dist, option, fig_pos)
% option: 1 = plot I, 2 = plot Theta

str_options = {'I(t)', 'Theta(t)'};

if isempty(cells_hist)
    msg = sprintf(' Unable to plot %s; ', str_options{option});
    return
end

N = size(cells_hist{1}, 1);
s = size(cells_hist{1}, 2);

tmax = numel(cells_hist)-1;
I = zeros(tmax+1, s);
Theta = zeros(tmax+1, s);

for i=1:tmax+1
    for j=1:s
        %cells = cells_hist{i}{j};
        cells = cells_hist{i};
        %[I(i,j), Theta(i,j)] = moranI(cells, a0*dist);
        [I(i,j), Theta(i,j)] = moranI(cells(:, j), a0*dist);
    end
end
%%        
h = figure; %(1+option);
cla(h, 'reset');
hold on
plot_clrs = [1 0 0;
             0 0 1];
if length(cells_hist) < 100
    lw = 1; ps = 'o-';
elseif length(cells_hist) < 500
    lw = 1; ps = '.-';
else
    lw = 0.5; ps = '.-';
end
          
for i=1:s
    switch option
        case 1
            y = I(:,i);
            h.Name = 'I_vs_t';
        case 2
            y = Theta(:,i);
            h.Name = 'Theta_vs_t';
    end
    clr = plot_clrs(i,:);
    plot(0:tmax, y, ps, 'LineWidth', lw, 'Color', clr);
end

%set(0, 'DefaultTextInterpreter', 'latex');
xlabel('$$t$$', 'Interpreter', 'latex');
plot_labels = {'$$I(t)$$', '$$\Theta(t)/f_N$$'};
ylabel(plot_labels{option}, 'Interpreter', 'latex');
legend(num2cell(string(1:s)));
set(gca, 'FontSize', 24);
xlim([0 tmax])
ylim([-1 1]);

set(h, 'Units', 'inches', 'Position', fig_pos);

msg = sprintf('Successfully plotted %s; ', str_options{option});  
end