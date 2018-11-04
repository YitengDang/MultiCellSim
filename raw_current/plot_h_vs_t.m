function [msg, h, h_all] = plot_h_vs_t(cells_hist, dist, a0, Con, K, noise, fig_pos)
% Only for 1 cell type

% msg : output message
% h : figure handle (not pseudo-energy!)
% h_all : pseudo-energy

h = [];
h_all = [];
if isempty(cells_hist)
    msg = 'Unable to plot h(t); ';
    return
end

s = size(cells_hist{1}, 2);
if s>1
    msg = 'Cannot plot h(t). More than 1 signalling molecule.';
    return
end

N = size(cells_hist{1}, 1);
tmax = numel(cells_hist)-1;
h_all = zeros(tmax+1, s);

% calculate interaction matrix
Rcell = 0.2*a0;
idx = dist>0;
M = ones(size(dist));
M(idx) = sinh(Rcell)./(a0*dist(idx)).*exp(Rcell-a0*dist(idx));

%
for i=1:tmax+1
    cells0 = cells_hist{i};
    for j=1:s
        cells = cells0(:, j);
        C0 = 1 + (Con-1).*cells;
        Y = M*C0;
        %dK = normrnd(0, noise, size(Y)); % Gaussian noise
        %K = K + dK;
        X = (2*cells-1);
        h_all(i,j) = -sum(X.*(Y-K));
    end
end
%%        
h = figure; %(6);
cla(h, 'reset');
hold on
plot_clrs = [1 0 0; 
                0 0 1];
if length(cells_hist) < 100
    ps = 'o-'; lw = 1;
elseif length(cells_hist) < 500
    ps = '.-'; lw = 1;
else
    ps = '.-'; lw = 0.5;
end
            
for i=1:s
    clr = plot_clrs(i,:);
    plot(0:tmax, h_all(:,i)/N, ps, 'LineWidth', lw, 'Color', clr);
end
xlabel('$$t$$', 'Interpreter', 'latex');
ylabel('$$h$$', 'Interpreter', 'latex');
legend(num2cell(string(1:s)));
set(gca, 'FontSize', 24);
xlim([0 tmax])
%ylim([0 1]);

set(h, 'Units', 'inches', 'Position', fig_pos);
h.Name = 'h_vs_t'; 

msg = 'Successfully plotted h(t); ';    
end