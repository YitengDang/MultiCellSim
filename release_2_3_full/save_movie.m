function save_movie(cells_hist, pos, cell_type, disp_mol, fname_out, frame_rate, rcell)
    %frames = struct('cdata',[],'colormap',[]);
    
    % Options
    %frame_rate = 5; % frames/second
    format = 'Motion JPEG AVI'; %movie format 
    % 'Motion JPEG AVI' <- default, works best
    % 'Uncompressed AVI' <- high quality(?), large file
    % 'MPEG-4' <- .mp4
    % 'Archival' <- unknown ext
    % 'Motion JPEG 2000' <- unknown ext
    
    % Save movie
    myVideo = VideoWriter(fname_out, format); %, 'Uncompressed AVI');
    myVideo.FrameRate = frame_rate;  % Default 30
    open(myVideo);
    
    % replay trajectory externally
    t = 0;
    h = figure(100);
    clf(h, 'reset');
    for t=1:length(cells_hist)
        cells = cells_hist{t};
        update_cell_figure_external(h, pos, cells, cell_type, t, disp_mol, rcell)                
        %frames(t) = getframe(h);
        frame = getframe(h);
        writeVideo(myVideo, frame);
    end
    
    cells = cells_hist{t};
    t = t+1;
    update_cell_figure_external(h, pos, cells, cell_type, t, disp_mol, rcell)                
    %frames(t) = getframe(h);
    frame = getframe(h);
    writeVideo(myVideo, frame);
        
    close(myVideo);
end