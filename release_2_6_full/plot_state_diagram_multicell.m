function [phase, A, h] = plot_state_diagram_multicell(gz, a0, rcell, M_int,...
    Con, Coff, K, lambda12, geneLogic)
% Plots the state diagram of the multicellular system with TWO types of
% signalling molecules

% Parameters
% lattice parameters
%gz = 15;
%N = gz^2;
%a0 = 1.5;
%rcell = 0.2;
Rcell = rcell*a0;

% circuit parameters
%Con = [18 16];
%Coff = [1 1];
%M_int = [1 1; -1 -1];
%K = [3 12; 13 20]; % K(i,j): sensitivity of type i to type j molecules
lambda = [1 lambda12]; % diffusion length (normalize first to 1)

% Exclude single cell
if gz==1
   [A, phase] = plot_state_diagram_onecell(M_int, Con, Coff, K, geneLogic);
   return 
end

% calculate fN
fN = zeros(2,1);
[dist, ~] = init_dist_hex(gz, gz);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence

fN(1) = sum(sinh(Rcell)*sum(exp((Rcell-r)./lambda(1)).*(lambda(1)./r)) ); % calculate signaling strength
fN(2) = sum(sinh(Rcell)*sum(exp((Rcell-r)./lambda(2)).*(lambda(2)./r)) ); % calculate signaling strength

if geneLogic == 1
    %% Calculate transition tables
    % Determine phases
    R1 = (repmat((1+fN').*Coff, 2, 1) - K) > 0; % Everything ON
    R2 = ((repmat(Con + fN'.*Coff, 2, 1) - K) > 0 & (repmat((1+fN').*Coff, 2, 1) - K) < 0); % ON remains ON & not all ON
    R3 = ((Coff + repmat(fN'.*Con, 2, 1) - K) < 0 & (repmat((1+fN').*Con, 2, 1) - K) > 0) ; % OFF remains OFF & not all OFF
    R4 = (repmat((1+fN').*Con, 2, 1) - K) < 0; % Everything OFF

    phase = 1 + R1 + 2*R2 + 3*R3 + 4*R4;
    % 1: none (activation-deactivation)
    % 2: all ON(+) / OFF(-)
    % 3: ON->ON/activation (+) / ON-> OFF (-) 
    % 4: OFF->OFF/deactivation (+) / OFF->ON (-)
    % 5: all OFF(+) / ON(-)
    % 6: autonomy (+) / autonomous oscillations (-)

    % Map from phase to diagram
    % state | activation/repression | input molecule (1/2)
    g_map = cell(2, 6, 2);
    % 0=OFF, 1:ON, 2:UNKNOWN
    % activation 
    g_map{1,1,1} = 2*ones(2);
    g_map{1,1,2} = 2*ones(2);
    g_map{1,2,1} = ones(2);
    g_map{1,2,2} = ones(2);
    g_map{1,3,1} = [2 2; 1 1];
    g_map{1,3,2} = [2 1; 2 1];
    g_map{1,4,1} = [0 0; 2 2];
    g_map{1,4,2} = [0 2; 0 2];
    g_map{1,5,1} = zeros(2);
    g_map{1,5,2} = zeros(2);
    g_map{1,6,1} = [0 0; 1 1];
    g_map{1,6,2} = [0 1; 0 1];
    % repression 
    %(note: this is precisely NOT g_map{1,:,:} in the three-val
    % boolean algebra with NOT 2 = 2)
    g_map{2,1,1} = 2*ones(2);
    g_map{2,1,2} = 2*ones(2);
    g_map{2,2,1} = zeros(2);
    g_map{2,2,2} = zeros(2);
    g_map{2,3,1} = [2 2; 0 0];
    g_map{2,3,2} = [2 0; 2 0];
    g_map{2,4,1} = [1 1; 2 2];
    g_map{2,4,2} = [1 2; 1 2];
    g_map{2,5,1} = ones(2);
    g_map{2,5,2} = ones(2);
    g_map{2,6,1} = [1 1; 0 0];
    g_map{2,6,2} = [1 0; 1 0];

    gij = cell(2);
    X_out = cell(2, 1);
    for i=1:2
        if all(M_int(i,:)==0)
            fprintf('No input for gene %d \n', i);
            % no input => output=initial state
            X1_in = [0 0; 1 1]; 
            X2_in = [0 1; 0 1];
            X_in = (i==1).*X1_in + (i==2).*X2_in; 
            X_out{i} = X_in;
        else
            % normal case
            for j=1:2
                if M_int(i,j)~=0
                    idx = (M_int(i,j)==1) + (M_int(i,j)==-1)*2;
                    gij{i,j} = g_map{idx, phase(i,j), j};
                else
                    gij{i,j} = ones(2); % Fixed ambiguous inputs
                end
            end
            X_out{i} = and3(gij{i,1}, gij{i,2}); % three-valued logic
        end
    end

    %% Display tables
    %{
    h1 = figure(1);

    % set colormap
    uniq_out = unique([X_out{1} X_out{2}]);
    colormap(map(uniq_out+1, :));

    subplot(1, 2, 1);
    imagesc([0 1], [0 1], X_out{1})
    set(gca, 'YDir', 'normal');
    xticks([0 1]);
    yticks([0 1]);
    set(gca, 'FontSize', 24);
    title('$$X^{(1)}_{out}$$ ')
    xlabel('$$X^{(1)}_{in}$$ ')
    ylabel('$$X^{(2)}_{in}$$ ')
    map = [1, 0, 0
        0, 1, 0
        0.9, 0.9, 0.1];
    cb = colorbar();
    caxis([min(uniq_out)-0.5 max(uniq_out)+0.5]);
    set(cb, 'YTick', uniq_out);

    subplot(1, 2, 2);
    imagesc([0 1], [0 1], X_out{2})
    set(gca, 'YDir', 'normal');
    xticks([0 1]);
    yticks([0 1]);
    set(gca, 'FontSize', 24);
    title('$$X^{(2)}_{out}$$ ')
    xlabel('$$X^{(1)}_{in}$$ ')
    ylabel('$$X^{(2)}_{in}$$ ')
    set(h1, 'Units', 'Inches', 'Position', [1 1 11 5]);
    cb = colorbar();
    caxis([min(uniq_out)-0.5 max(uniq_out)+0.5]);
    set(cb, 'YTick', uniq_out);
    %}
    %% Calculate state diagram
    A = zeros(4); % graph adjacency matrix
    for i=1:2
        for j=1:2
            state_in = i + 2*(j-1); 
            X_out_this = [X_out{1}(i,j) X_out{2}(i,j)]; % tentative 
            %disp(X_out_this)
            if all(X_out_this~=2) % unambiguous out state
                state_out = X_out_this(1)+1 + 2*X_out_this(2); % (i,j) -> idx
                A(state_in, state_out) = 1;
                %disp(state_out);
            elseif sum(X_out_this==2)==1 % semi-definite
                if (X_out_this(1)==2)
                    X_out_both = [0 X_out_this(2); 1 X_out_this(2)];
                elseif (X_out_this(2)==2)
                    X_out_both = [X_out_this(1) 0; X_out_this(1) 1];
                end
                state_out = X_out_both*[1; 2]+1;
                %[X_out_both(1,1)+1 + 2*X_out_both(1,2);...
                %    X_out_both(2,1)+1 + 2*X_out_both(2,2)];
                A(state_in, state_out) = 1;
                %disp(state_out);
            elseif sum(X_out_this==2)==2 
                A(state_in, :) = 1;
            end
        end
    end
elseif geneLogic == 2 % OR Logic
    disp('OR Gene Logic: not yet implemented');
    phase = []; 
    A = [];
    h = [];
    return
else
    error('Invalid Gene Logic (choose AND or OR logic)');
end

%% Draw state diagram
h = figure; %(10);
h.Name = 'plot_state_diagram_multicell';
hold on
s = [0 1 0 1];
t = [0 0 1 1];
%A = ones(4);
Gs = digraph(A);
nLabels = {};
g=plot(Gs, 'XData', s, 'YData', t, 'ArrowSize', 20, 'EdgeAlpha', 1, ...
    'LineWidth', 3, 'EdgeColor', 'k',...
    'Marker', 'o', 'MarkerSize', 100, 'NodeColor', [0.2 0.2 0.2], 'NodeLabel', nLabels);
% Make edges dashed if state has two outgoing edges
for i=1:4
    if sum(A(i,:))==2
        idx = find(A(i,:));
        highlight(g, i, idx, 'LineStyle', '--', 'LineWidth', 2);
    elseif sum(A(i,:))==4
        highlight(g, i, 1:4, 'LineStyle', ':', 'LineWidth', 2);
    end
end
text(s-0.11,t+0.019,{'(0,0)','(1,0)','(0,1)','(1,1)'}, 'Color', 'w', 'FontSize', 32)
ax = gca;
axis([-0.4 1.4 -0.4 1.4]);
ax.Visible = 'off';
h.Color = [1 1 1];
%set(ax, 'Units', 'Inches', 'Position', [0 0 9 8]);
%set(h2, 'Units', 'Inches', 'Position', [1 1 9 8]);
set(ax, 'Units', 'Inches', 'Position', [0 0 7 6]);
set(h, 'Units', 'Inches', 'Position', [0.2 0.2 7 6]);



end

%% Functions
% 3-valued AND function (see three_valued_logic.m)
function out = and3(x,y)
    out = min(x.*y, 2);
end