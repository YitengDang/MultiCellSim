clear all
close all
%% Load file
folder = fullfile(pwd, 'data', 'time_evolution');
[file, path] = uigetfile(fullfile(folder, '\*.mat'), 'Load saved simulation');

if isnumeric(file) && isnumeric(path) % if user clicked "cancel"
    %app.MessagesTextArea.Value = sprintf('Simulation not loaded.');
    return
end
load(fullfile(path, file));
%% Plot figure
fig_pos = [1 1 8 6];
t1 = 40;
cells_hist_plot = cells_hist(t1+1:end);
msg = plot_p_vs_t(cells_hist_plot, fig_pos)

%% Save new figure
folder = fullfile(pwd, 'data', 'time_evolution', 'sample_trajectories');
ext = '.pdf';
num_fig = 6; % total number of figures (needs to be updated)
            
fig_label = cell(num_fig, 1);
fig_label{1} = 'p_vs_t';
fig_label{2} = 'I_vs_t';
fig_label{3} = 'Theta_vs_t';
fig_label{4} = 'p_I_vs_t';
fig_label{5} = 'p_Theta_vs_t';
fig_label{6} = 'h_vs_t';

for i=1:num_fig
    h = figure(i);
    if ~isempty(h.Children)
        fname = fullfile(folder, strcat('temp', '_', fig_label{i}, ext));
        [file, path] = uiputfile(fname);

        if ~file % if user pressed cancel, continue to next fig
            continue;
        end
        [~, name, ext] = fileparts(file); %separate file name from extension

        fname_out = fullfile(path, name);
        save_figure(h, 10, 8, fname_out, ext);
    else
        close(h);
    end
end