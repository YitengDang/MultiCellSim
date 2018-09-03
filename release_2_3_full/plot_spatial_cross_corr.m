%{
clear all
close all
set(0, 'defaulttextinterpreter', 'latex');
folder = 'D:\Multicellularity\data\two_signals\time_evolution\temp';
fname_str = 'test3';
load(fullfile(folder, fname_str));
%%
cells = cells_hist{1};
N=size(cells, 1);
dH = sum(cells(:,1) ~= cells(:,2))/N;

idx = find(cells(:,1) | cells(:,2));
J = sum(cells(idx,1) & cells(idx,2))/numel(idx);
Jd = 1-J;

[dH_all, dJ_all] = plot_spatial_cross_corr_temp(cells_hist);
%}
function [dH, dJ] = plot_spatial_cross_corr(cells_hist, fig_pos)
    if nargin<2
        fig_pos = [0.1 0.1 7 5];
    end
    set(0, 'defaulttextinterpreter', 'latex');
    N = size(cells_hist{1}, 1);
    n_dat = numel(cells_hist);
    
    dH = zeros(n_dat, 1);
    dJ = zeros(n_dat, 1);
    for i=1:n_dat
        cells = cells_hist{i};
        
        % (1) Hamming dist.
        % NB same as counting (0,1) and (1,0) together
        dH(i) = sum(cells(:,1) ~= cells(:,2))/N;

        % (2) Jaccard index/coeff = 1 - (Jacard distance)
        idx = find(cells(:,1) | cells(:,2));
        if ~isempty(idx)
            J = sum(cells(idx,1) & cells(idx,2))/numel(idx);
            dJ(i) = 1-J;
        else
            dJ(i) = 0;
        end
        
        
    end
    
    %% Plot result
    h=figure;
    hold on
    plot(0:n_dat-1, dH, 'LineWidth', 1.5);
    plot(0:n_dat-1, dJ, 'LineWidth', 1.5);
    legend({'Hamming', 'Jaccard'});
    xlabel('time');
    ylabel('distance');
    set(gca, 'FontSize', 24);
    set(h, 'Units', 'inches', 'Position', fig_pos);
    xlim([0 n_dat-1]);
    ylim([0 1]);
    % correlation between Hamming and Jacard distance?
    %{
    figure;
    scatter(dH, dJ);
    [corr_dH_dJ, pval] = corr(dH, dJ);
    title(sprintf('$$\\rho = %.3f, p = %.3f$$', corr_dH_dJ, pval));
    xlabel('Hamming dist');
    ylabel('Jaccard dist');
    %}
end