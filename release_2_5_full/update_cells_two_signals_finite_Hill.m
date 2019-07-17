function [cells_out, changed] = ...
    update_cells_two_signals_finite_Hill(cells, dist, M_int, a0,...
    Rcell, Con, Coff, K, lambda, hill, noise, logic)
% Update cells without noise in a positive feedback loop with infinite hill
% coefficient
% noise: fractional noise term, represents the width of the normal
% distribution as a fraction of the mean of the noisy parameter
% logic: specifies logic gate for integrating multiple signals. 
% Options: AND (1, default), OR (2).
if nargin<12
    logic = 1; % AND
else
    if logic~=1 && logic~=2
        error('Invalid logic function specified');
    end
end
%%
if noise>1
    error('Invalid noise term. Input: %.2f, but 0 <= noise <= 1 is required', noise);
end
N = size(cells, 1);

% Account for self-influence
idx = dist>0;
% TO DO: vectorize / combine
M1 = ones(size(dist)); 
M1(idx) = sinh(Rcell)./(a0*dist(idx)/lambda(1))...
    .*exp((Rcell-a0*dist(idx))/lambda(1));
M2 = ones(size(dist)); 
M2(idx) = sinh(Rcell)./(a0*dist(idx)/lambda(2))...
    .*exp((Rcell-a0*dist(idx))/lambda(2));

% Concentration in each cell
C0 = Coff + (Con-Coff).*cells; 

% Reading of each cell
Y1 = M1*C0(:, 1); 
Y2 = M2*C0(:, 2);
Y = [Y1 Y2];

% Add noise to K
K_cells = K.*ones(2, 2, N);
%if noise>0
dK = K_cells.*normrnd(0, noise, 2, 2, N); 
K_cells = max(K_cells + dK, 1); % do not allow K < 1 = Coff
%end
%%
% Multiplicative interaction
if hill==Inf
    out = zeros(N, 2,2);
    out(:,1,1) = ((Y1-squeeze(K_cells(1,1,:)))*M_int(1,1) > 0);
    out(:,1,2) = ((Y2-squeeze(K_cells(1,2,:)))*M_int(1,2) > 0);
    out(:,2,1) = ((Y1-squeeze(K_cells(2,1,:)))*M_int(2,1) > 0);
    out(:,2,2) = ((Y2-squeeze(K_cells(2,2,:)))*M_int(2,2) > 0);
    
    switch logic
        case 1 % AND logic
            out = out + (1 - permute(repmat(abs(M_int), 1, 1, N), [3 2 1]) );
            X1 = out(:,1,1) & out(:,1,2);
            X2 = out(:,2,1) & out(:,2,2);
        case 2 % OR logic
            X1 = out(:,1,1) | out(:,1,2);
            X2 = out(:,2,1) | out(:,2,2);
    end
elseif hill > 0
    K_cells_1 = squeeze(K_cells(1,:,:))';
    K_cells_2 = squeeze(K_cells(2,:,:))';
    
    fX1 = 1./( 1 + ...
        ((Y./K_cells_1).*(1-M_int(1,:))/2).^hill + ... % case repression 
        ((K_cells_1./Y).*(1+M_int(1,:))/2).^hill ... % case activation
        ).*abs(M_int(1,:)) + ... 
        (1-abs(M_int(1,:))).*ones(N,2); % case no interaction
    fX2 = 1./( 1 + ...
        ((Y./K_cells_2).*(1-M_int(2,:))/2).^hill + ... % case repression 
        ((K_cells_2./Y).*(1+M_int(2,:))/2).^hill ... % case activation
        ).*abs(M_int(2,:)) + ... 
        (1-abs(M_int(2,:))).*ones(N,2); % case no interaction
    
    switch logic
        case 1 % AND logic
            X1 = prod(fX1, 2);
            X2 = prod(fX2, 2);
        case 2 % OR logic 
            % (A OR B) = NOT( NOT(A) AND NOT(B) )
            X1 = 1 - (1- fX1(:, 1)).*(1- fX1(:, 2));
            X2 = 1 - (1- fX2(:, 1)).*(1- fX2(:, 2));
    end
    
end

%% output
cells_out = [X1 X2];
% if no connections to a gene, output = input (remains constant)
idx2 = find(sum(abs(M_int), 2)==0); % find channel(s) that don't have any input
cells_out(:, idx2) = cells(:, idx2); % revert to input

changed = ~isequal(cells_out, cells);

end