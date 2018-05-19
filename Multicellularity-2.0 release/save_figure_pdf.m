function save_figure_pdf(h_fig, width, height, path_out)
% Adjust the height and width of a figure and save it in PDF
    set(h_fig,'Units','Inches');
    set(h_fig, 'Position', [0 0 width height ])
    pos = get(h_fig,'Position');
    set(h_fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h_fig, path_out,'-dpdf','-r0')
end