function [cells_out, changed] = ...
    update_cells_two_signals(cells, dist, M_int, a0, Rcell, Con, K, lambda, noise)
% Update cells without noise in a positive feedback loop with infinite hill
% coefficient
% Multiplicative interaction

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
C0 = 1 + (Con-1).*cells; 

% Reading of each cell
Y1 = M1*C0(:, 1); 
Y2 = M2*C0(:, 2);

% Add noise to K
dK = normrnd(0, noise, 2, 2);
K = K + dK;

% Multiplicative interaction (AND)
out11 = ((Y1-K(1,1))*M_int(1,1) > 0) + (1 - abs(M_int(1,1)));
out12 = ((Y2-K(1,2))*M_int(1,2) > 0) + (1 - abs(M_int(1,2)));
out21 = ((Y1-K(2,1))*M_int(2,1) > 0) + (1 - abs(M_int(2,1)));
out22 = ((Y2-K(2,2))*M_int(2,2) > 0) + (1 - abs(M_int(2,2)));

cells_out = [out11.*out12 out21.*out22];

% Additive interaction (OR)


changed = ~isequal(cells_out, cells);
%%
% hamiltonian
%hlist = -(2*cells-1).*([Y1 Y2]-K);
%h = sum(hlist);        