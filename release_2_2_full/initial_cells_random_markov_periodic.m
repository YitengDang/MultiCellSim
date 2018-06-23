function [pos, dist, fN0, rejections] = initial_cells_random_markov_periodic(n, Lx, R, mcsteps, nodisplay)
% Note: even perfect arrangement might not be accepted because distances
% are rounded off.
% Places cells randomly in a continuous space
%{
clear all
L = 1;
R = 0.02;
n = round(L/R/5); % nmax = L/R
%}
if nargin<5
    nodisplay = 0;
end

if ~nodisplay
    disp('Initiating initial lattice...');
end
%---------main code below---------------
N = n^2;

% hexagonal placement
delx = Lx/n;
dely = sqrt(3)/2*delx;
Ly = dely*n;
[xm, ym] = meshgrid(0:n-1, 0:n-1);
x = (xm+mod(ym,2)/2)*delx;
y = ym*dely;

% store initial config
%x_ini = x;
%y_ini = y;

% Calculate distance
dist = calc_dist_periodic(x, y, Lx, Ly);
dist = round(dist, 10);
%pos = [x(:) y(:)];

% Calculate default fN
fN = zeros(N, 1);
for i=1:N
    dist_vec = dist(i,:);
    r = dist_vec(dist_vec>0); % exclude self influence
    fN(i) = sum(sinh(R)*sum(exp(R-r)./r)); % calculate signaling strength
end
fN0 = mean(fN);
if ~nodisplay
    fprintf('fN0 (default fN) = %.2f \n', fN0);
end
%%
% Key parameters
% mcsteps = 10^5; % Monte Carlo steps
delta = (delx - 2*R)/4; % step size
step = 0;
rejections = 0;
while step < mcsteps
    %disp(step)
    if mod(step, 10^4)==0 && ~nodisplay
        fprintf('MC steps: %d, ', step);
        fprintf('Rejections: %d \n', rejections);
    end
    % Random step
    j = randi(N); %selected particle
    
    x_new = x; 
    y_new = y;
    x_new(j) = mod(x_new(j) + delta*rand(), Lx);
    y_new(j) = mod(y_new(j) + delta*rand(), Ly);
    
    % Calculate distance
    dist_new = calc_dist_periodic(x_new, y_new, Lx, Ly);

    % Check distance
    if all(dist_new(dist_new>0) >= 2*R) 
        x = x_new;
        y = y_new;
        dist = round(dist_new, 10);
        step = step+1;
    else
        rejections = rejections + 1;
        %fprintf('Rejected! \n');
    end
end

pos = [x(:) y(:)];
% adjust distances to scale so that nearest neighbour distances are
% 1 (in case of a perfect lattice); without normalization this
% would be Lx/n.
dist = dist/(Lx/n); 

if ~nodisplay
    fprintf('final MC trials: %d \n', mcsteps);
    fprintf('final rejections: %d \n', rejections);
end
%-----------end main code----------------------
%% Draw configuration
%{
hin = figure();
clf(hin,'reset');
title(sprintf('N = %d, R = %.2f', N, R), ...
    'FontSize', 24);
%set(gca,'YTick',[],'XTick',[]);
set(gca,'DataAspectRatio', [1 1 1]);
axis([0 L 0 L]);
box on
hold on
for i=1:N
    position = [x(i)-R y(i)-R 2*R 2*R];
    face_clr = 'k';
    curv = [1 1];
    rectangle('Position', position, 'FaceColor', face_clr, ...
        'EdgeColor', 'k', 'Curvature', curv);
    % draw another circle for pbc
    cond = [x(i)<R 1-x(i)<R y(i)<R 1-y(i)<R];
        % sides
    if cond(1)
        rectangle('Position', [L+x(i)-R y(i)-R 2*R 2*R], 'FaceColor', face_clr, ...
        'EdgeColor', face_clr, 'Curvature', curv);
    elseif cond(2)
        rectangle('Position', [-L+(x(i)-R) y(i)-R 2*R 2*R], 'FaceColor', face_clr, ...
        'EdgeColor', face_clr, 'Curvature', curv);
    end
    if cond(3)
        rectangle('Position', [x(i)-R L+y(i)-R 2*R 2*R], 'FaceColor', face_clr, ...
        'EdgeColor', face_clr, 'Curvature', curv);
    elseif cond(4)
        rectangle('Position', [x(i)-R -L+(y(i)-R) 2*R 2*R], 'FaceColor', face_clr, ...
        'EdgeColor', face_clr, 'Curvature', curv);
    end
    % corners
    if cond(1) && cond(3)
        rectangle('Position', [L+x(i)-R L+(y(i)-R) 2*R 2*R], 'FaceColor', face_clr, ...
        'EdgeColor', face_clr, 'Curvature', curv);
    elseif cond(2) && cond(4)
        rectangle('Position', [-L+(x(i)-R) -L+(y(i)-R) 2*R 2*R], 'FaceColor', face_clr, ...
        'EdgeColor', face_clr, 'Curvature', curv);
    elseif cond(1) && cond(4)
        rectangle('Position', [L+(x(i)-R) -L+(y(i)-R) 2*R 2*R], 'FaceColor', face_clr, ...
        'EdgeColor', face_clr, 'Curvature', curv);
    elseif cond(2) && cond(3)
        rectangle('Position', [-L+(x(i)-R) L+(y(i)-R) 2*R 2*R], 'FaceColor', face_clr, ...
        'EdgeColor', face_clr, 'Curvature', curv);    
    end
end
hold off
drawnow;
%}
%% Check whether final configuration is different enough from initial config
% Delta r = distance between initial and final position of particle
%{
dr = sqrt( (x-x_ini).^2  + (y-y_ini).^2 );
h=figure();
histogram(dr/L, 'Normalization', 'pdf');
title( strcat('$$\langle \Delta r / L \rangle$$ = ', sprintf('%.2f', mean(mean(dr)))), 'Interpreter', 'latex');
xlabel('\Delta r / L');
ylabel('Probability');
set(gca, 'FontSize', 16);
set(h, 'Position', [500 200 840/1.5 720/1.5]);
%}
%% Test whether distribution is random enough
%{
dist2 = reshape(dist(dist>0), [N-1 N]);
nnd = min(dist2);
nnd_av = mean(nnd);

h=figure();
hold on
histogram(nnd, 'Normalization', 'pdf');
xlabel('NND');
ylabel('Probability');
title(sprintf('%s NND %s = %.2f', '\langle', '\rangle', nnd_av));
set(gca, 'FontSize', 24);
set(h, 'Position', [500 200 840 720]);

% Clark & Evans calculation
rho = @(r) 2*pi*N/L^2*r.*exp(-pi*N/L^2*r.^2);
rvals = linspace(0, max(nnd)*1.2, 1000);
plot(rvals, rho(rvals), 'LineWidth', 2);
%xlim([0 1]);

% set image properties
%h = gcf;
%set(h,'Units','px');
%set(h, 'Position', [500 300 600 600]);
%}