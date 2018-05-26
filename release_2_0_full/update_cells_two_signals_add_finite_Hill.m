function [cells_out, changed] = ...
    update_cells_two_signals_add_finite_Hill(cells, dist, M_int, a0, Rcell,...
    Con, Coff, K, lambda, hill, noise)
% Update cells without noise in a positive feedback loop with infinite hill
% coefficient
N = size(cells, 1);

% correct for non-zero entries in non-existing interactions
Coff = Coff.*abs(M_int);
Con = Con.*abs(M_int);
% sum production rates due to different molecules
Coff_tot = sum(Coff, 2); 
Con_tot = sum(Con, 2);

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
C0 = sum(Con-Coff, 2)'.*cells + sum(Coff,2)';

% Reading of each cell
Y1 = M1*C0(:, 1); % molecule 1
Y2 = M2*C0(:, 2); % molecule 2
Y = [Y1 Y2];
%%
% Add noise to K
K_cells = K.*ones(2, 2, N);
dK = normrnd(0, noise, 2, 2, N);
K_cells = max(K_cells + dK, 1); % do not allow K < 1 = Coff
%%
% explicit calculation
%S11 = (Con(1,1)-Coff(1,1)).*fX1(:,1) + Coff(1,1);
%S12 = (Con(1,2)-Coff(1,2)).*fX1(:,2) + Coff(1,2);
%S1b = S11 + S12;

% vectorized calculation 
% factors 1+M_int/2 and 1-M_int/2 to distinguish activation and repression
if hill==Inf
    %fX1 = ((Y > K(1,:)).*(1+M_int(1,:))/2 + (Y < K(1,:)).*(1-M_int(1,:))/2).*abs(M_int(1,:));
    %fX2 = ((Y > K(2,:)).*(1+M_int(2,:))/2 + (Y < K(2,:)).*(1-M_int(2,:))/2).*abs(M_int(2,:));
    fX1 = ((Y > squeeze(K_cells(1,:,:))').*(1+M_int(1,:))/2 + (Y < squeeze(K_cells(1,:,:))').*(1-M_int(1,:))/2).*abs(M_int(1,:));
    fX2 = ((Y > squeeze(K_cells(2,:,:))').*(1+M_int(2,:))/2 + (Y < squeeze(K_cells(2,:,:))').*(1-M_int(2,:))/2).*abs(M_int(2,:));
elseif hill>0
    %fX1 = (Y.^hill.*(1+M_int(1,:))/2 + (K(1,:).^hill.*(1-M_int(1,:))/2))...
    %    ./(K(1,:).^hill+Y.^hill).*abs(M_int(1,:));
    %fX2 = (Y.^hill.*(1+M_int(2,:))/2 + (K(2,:).^hill.*(1-M_int(2,:))/2))...
    %    ./(K(2,:).^hill+Y.^hill).*abs(M_int(2,:));
    fX1 = (Y.^hill.*(1+M_int(1,:))/2 + (squeeze(K_cells(1,:,:))'.^hill.*(1-M_int(1,:))/2))...
        ./(squeeze(K_cells(1,:,:))'.^hill+Y.^hill).*abs(M_int(1,:));
    fX2 = (Y.^hill.*(1+M_int(2,:))/2 + (squeeze(K_cells(2,:,:))'.^hill.*(1-M_int(2,:))/2))...
        ./(squeeze(K_cells(2,:,:))'.^hill+Y.^hill).*abs(M_int(2,:));
else 
    warning('Hill coefficient ill-defined');
    return
end

S1 = sum(((Con(1,:)-Coff(1,:)).*fX1 + Coff(1,:)), 2); % vectorized
S2 = sum(((Con(2,:)-Coff(2,:)).*fX2 + Coff(2,:)), 2); % vectorized

X1 = (S1 - Coff_tot(1))./(sum(Con_tot(1), 2) - sum(Coff_tot(1), 2));
X2 = (S2 - Coff_tot(2))./(sum(Con_tot(2), 2) - sum(Coff_tot(2), 2));

%% output
cells_out = [X1 X2];
changed = ~isequal(cells_out, cells);
%%
% hamiltonian
%hlist = -(2*cells-1).*([Y1 Y2]-K);
%h = sum(hlist);        