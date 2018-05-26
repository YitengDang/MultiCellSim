function [cells_out, changed] = ...
    update_cells_two_signals_multiply_finite_Hill(cells, dist, M_int, a0,...
    Rcell, Con, Coff, K, lambda, hill, noise)
% Update cells without noise in a positive feedback loop with infinite hill
% coefficient

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
dK = normrnd(0, noise, 2, 2);
K = K + dK;
%%
% Multiplicative interaction
if hill==Inf
    out11 = ((Y1-K(1,1))*M_int(1,1) > 0) + (1 - abs(M_int(1,1)));
    out12 = ((Y2-K(1,2))*M_int(1,2) > 0) + (1 - abs(M_int(1,2)));
    out21 = ((Y1-K(2,1))*M_int(2,1) > 0) + (1 - abs(M_int(2,1)));
    out22 = ((Y2-K(2,2))*M_int(2,2) > 0) + (1 - abs(M_int(2,2)));
    X1 = out11.*out12;
    X2 = out21.*out22;
elseif hill > 0
    fX1 = (Y.^hill.*(1+M_int(1,:))/2 + (K(1,:).^hill.*(1-M_int(1,:))/2))...
        ./(K(1,:).^hill+Y.^hill).*abs(M_int(1,:)) + (1-abs(M_int(1,:))).*ones(N,2);
    fX2 = (Y.^hill.*(1+M_int(2,:))/2 + (K(2,:).^hill.*(1-M_int(2,:))/2))...
        ./(K(2,:).^hill+Y.^hill).*abs(M_int(2,:)) + (1-abs(M_int(2,:))).*ones(N,2);
    X1 = prod(fX1, 2);
    X2 = prod(fX2, 2);
end

%% output
cells_out = [X1 X2];
changed = ~isequal(cells_out, cells);
%%
% hamiltonian
%hlist = -(2*cells-1).*([Y1 Y2]-K);
%h = sum(hlist);        