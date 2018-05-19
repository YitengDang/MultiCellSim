function msg = plot_p_I_vs_t(cells_hist, dist, a0, Con, K, fN, option, fig_pos)
% option: 1 = plot (p,I), 2 = plot (p,Theta)

str_options = {'(p(t), I(t))', '(p(t), Theta(t))'};

if isempty(cells_hist)
    msg = sprintf(' Unable to plot %s; ', str_options{option});
    return
end

s = size(cells_hist{1}, 2);
N = size(cells_hist{1}, 1);
%if s>1
%    N = size(cells_hist{1}{1}, 1);
%else
%    N = size(cells_hist{1}, 1);
%end

tmax = numel(cells_hist)-1;
p = zeros(tmax+1, s);
I = zeros(tmax+1, s);
Theta = zeros(tmax+1, s);

for i=1:tmax+1
    if s>1
        for j=1:s
            %cells = cells_hist{i}{j};
            cells = cells_hist{i}; cells = cells(:, j);
            p(i,j) = sum(cells)/N;
            [I(i,j), Theta(i,j)] = moranI(cells, a0*dist);
        end
    else
        cells = cells_hist{i};
        p(i,1) = sum(cells)/N;
        [I(i,1), Theta(i,1)] = moranI(cells, a0*dist);
    end
end
%%
h4 = figure(option+3);
cla(h4, 'reset');
hold on

% plot hamiltonian first
if s==1
    pv = (0:N)/N;
    switch option
        case 1
            E = @(p, I) -0.5*(Con-1)*(1 + 4*fN.*p.*(1-p).*I + fN*(2*p-1).^2) ...
                -(2*p-1).*(0.5*(Con+1)*(1+fN) - K);
            yv = -0.15:0.05:1; % y = I
        case 2
            E = @(p, Theta) -0.5*(Con-1)*(1 + fN*Theta)...
                -(2*p-1).*(0.5*(Con+1)*(1+fN) - K);
            yv = -1:0.05:1; % y = Theta
    end
    [p_i,ymesh] = meshgrid(pv, yv);
    contourf(pv, yv, E(p_i, ymesh),'LineStyle', 'none')
    colormap(h4, 'summer')
    c = colorbar();
    c.Label.String = 'h'; 
    c.Label.FontSize = 24;
else
    ylim_vals = [-0.15 1]; %limits of ylim
    yv = ylim_vals;
end

% plot trajectory / trajectories
if length(cells_hist) < 100
    lw = 1.5;
elseif length(cells_hist) < 500
    lw = 1; 
else
    lw = 0.5; 
end

%plot_colours = get(gca,'colororder');
plot_clrs = [1 0 0; 
             0 0 1];
plots = zeros(s, 1);

% what to plot
switch option
    case 1
        y = I;
    case 2
        y = Theta;
        % plot Theta/fN = (2p-1)^2 (solution for uniform lattice
        pv = (0:N)/N;
        plot(pv, (2.*pv-1).^2, 'k--', 'LineWidth', 0.8)
end

% plot for all signals/genes
for i=1:s
    clr = plot_clrs(i,:);
    plot1 = plot(p(:,i), y(:,i), '--', 'LineWidth', lw, 'Color', clr);
    
    % store for legend
    plots(i) = plot1;
end

% plot start and end points
for i=1:s
    clr = plot_clrs(i,:);
    plot(p(1,i), y(1,i), 'o', 'MarkerSize', 8, 'Color', clr, 'MarkerFaceColor', clr);
    plot(p(end,i), y(end,i), 'o', 'MarkerSize', 8, 'Color', clr, 'MarkerFaceColor', 'w');
    plot(p(end,i), y(end,i), 'x', 'MarkerSize', 8, 'Color', clr);
end

%set(0, 'DefaultTextInterpreter', 'latex');
xlabel('$$p(t)$$', 'Interpreter', 'latex');
plot_labels = {'$$I(t)$$', '$$\Theta(t)/f_N$$'};
ylabel(plot_labels{option}, 'Interpreter', 'latex');
legend(plots, num2cell(string(1:s)) );
set(gca, 'FontSize', 24);
xlim([0 1])

ylim([yv(1) yv(end)]);

% set position
set(h4, 'Units', 'inches', 'Position', fig_pos);

msg = sprintf('Successfully plotted %s; ', str_options{option});  
end