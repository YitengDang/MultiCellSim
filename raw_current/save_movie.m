function save_movie(cells_hist, rcell, pos, pos_hist, disp_mol, fname_out,...
    frame_rate, t_ini, t_out)
    %frames = struct('cdata',[],'colormap',[]);
    
    if nargin<8
        t_ini = 0;
        t_out = length(cells_hist)-1;
    end
    % Options
    %frame_rate = 5; % frames/second
    format = 'Motion JPEG AVI'; %movie format 
    % 'Motion JPEG AVI' <- default, works best
    % 'Uncompressed AVI' <- high quality(?), large file
    % 'MPEG-4' <- .mp4
    % 'Archival' <- unknown ext
    % 'Motion JPEG 2000' <- unknown ext
    
    % Initiate movie
    myVideo = VideoWriter(fname_out, format); %, 'Uncompressed AVI');
    myVideo.FrameRate = frame_rate;  % Default 30
    open(myVideo);
    
    % replay trajectory externally
    %tt = 0;
    h = figure;
    clf(h, 'reset');
    
    if isempty(pos_hist)
        % same positions every time step
        [h_cells, h_borders] = reset_cell_figure(h, pos, rcell);
        for tt=t_ini+1:t_out+1
            cells = cells_hist{tt};
            update_cell_figure_external(h_cells, h_borders, cells, tt-1, disp_mol, pos);
            %update_cell_figure_external(h, pos, cells, cell_type, tt-1, disp_mol, rcell)                
            %frames(t) = getframe(h);
            frame = getframe(h);
            writeVideo(myVideo, frame);
        end
    else
        % new positions every time step
        [h_cells, h_borders] = reset_cell_figure(h, pos, rcell);
        for tt=t_ini+1:t_out+1
            cells = cells_hist{tt};
            pos = pos_hist{tt};
            update_cell_figure_external(h_cells, h_borders, cells, tt-1, disp_mol, pos);
            %update_cell_figure_external(h, pos, cells, cell_type, tt-1, disp_mol, rcell)                
            %frames(t) = getframe(h);
            frame = getframe(h);
            writeVideo(myVideo, frame);
        end
    end
    
    % add final state to show equilibrium (not applicable if not in
    % equilibrium)
    %{
    cells = cells_hist{t};
    t = t+1;
    update_cell_figure_external(h, pos, cells, cell_type, t-1, disp_mol, rcell)     
    frames(t) = getframe(h);
    frame = getframe(h);
    writeVideo(myVideo, frame);
    %}
        
    close(myVideo);
end