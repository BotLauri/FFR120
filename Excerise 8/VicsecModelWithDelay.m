%% Vicsec model with delay.
clear all
tic
% Setup.
N = 100; % Number of particles. 
L = 1000; % Side length. 
Rf = 20; % Interaction radius. 
eta = 0.4; % Noise level.
deltaT = 1; % Time step.
v = 3; % Speed.
rho = N/L^2; % Particle density.
S = 10^4; % Time steps.

%r = InitializePositions(N, L); % Random initialization. 
%save('initialRDelay', 'r');
load('initialRDelay', 'r');
initialR = r;
load('initialTheta', 'theta');

% Exercise 8.8
%H = [0:25];
H = [-1:0];
globAlignCoeff = zeros(length(H), S);
globClustCoeff = zeros(length(H), S);
for h = 1:length(H)
    
oldThetas = {}; oldThetas{end + 1} = {theta};
for m = 1:S
    isNeighbour = zeros(N) ~= 0;
    for i = 1:N
        isNeighbour(i, :) = FindNeighbours(r, i, Rf, L);
    end
    
    velocities = UpdateVelocities(v, theta);
    deltaR = velocities.*deltaT;
    r = UpdatePositions(r, deltaR, L);
    theta = UpdateOrientationWithDelay(theta, isNeighbour, eta, deltaT, H(h), oldThetas);
    oldThetas{end + 1} = {theta};
    globAlignCoeff(m) = CalculateGlobalAlignmentCoefficient(velocities, v);
    globClustCoeff(m) = CalculateGlobalClusteringCoefficent(r, Rf, L);
end

if (H(h) == 0)
    rh0 = r;
end

if (H(h) == 2)
    rh2 = r;
end

if (H(h) == 5)
    rh5 = r;
end

if (H(h) == 10)
    rh10 = r;
end

if (H(h) == 25)
    rh25 = r;
end

disp(H(h))

end

subplot(2, 3, 1)
PlotVoronoiDiagram(initialR, L)
axis equal
axis([-3/2*L 3/2*L -3/2*L 3/2*L])
title('Initial distribution.')

subplot(2, 3, 2)
PlotVoronoiDiagram(rh0, L)
axis equal
axis([-3/2*L 3/2*L -3/2*L 3/2*L])
title('Distribution after 10^4 timesteps for h = 0.')

subplot(2, 3, 3)
PlotVoronoiDiagram(rh2, L)
axis equal
axis([-3/2*L 3/2*L -3/2*L 3/2*L])
title('Distribution after 10^4 timesteps for h = 2.')

subplot(2, 3, 4)
PlotVoronoiDiagram(rh5, L)
axis equal
axis([-3/2*L 3/2*L -3/2*L 3/2*L])
title('Distribution after 10^4 timesteps for h = 5.')

subplot(2, 3, 5)
PlotVoronoiDiagram(rh10, L)
axis equal
axis([-3/2*L 3/2*L -3/2*L 3/2*L])
title('Distribution after 10^4 timesteps for h = 10.')

subplot(2, 3, 6)
PlotVoronoiDiagram(rh25, L)
axis equal
axis([-3/2*L 3/2*L -3/2*L 3/2*L])
title('Distribution after 10^4 timesteps for h = 25.')

%save('globalAlignmentCoefficient', 'globAlignCoeff');
%save('globalClusteringCoefficient', 'globClustCoeff');

toc
