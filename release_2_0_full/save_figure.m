function save_figure(h_fig, width, height, path_out, ext)
% Adjust the height and width of a figure and save it in PDF
    set(h_fig,'Units','Inches');
    set(h_fig, 'Position', [0 0 width height ])
    pos = get(h_fig,'Position');
    set(h_fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    if strcmp(ext, '.pdf')
        print(h_fig, path_out,'-dpdf','-r0');
    elseif strcmp(ext, '.eps')
        print(h_fig, path_out,'-depsc','-r0');
    elseif strcmp(ext, '.svg')
        print(h_fig, path_out,'-dsvg','-r0');
    end
end