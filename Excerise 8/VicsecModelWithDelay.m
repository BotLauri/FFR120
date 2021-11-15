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
S = 5*10^3; % Time steps.

r = InitializePositions(N, L); % Random initialization. 
save('initialRDelay', 'r');
load('initialRDelay', 'r');
initialR = r;
load('initialTheta', 'theta');
globAlignCoeff = zeros(1, S);
globClustCoeff = zeros(1, S);
isNeighbour = zeros(N) ~= 0;
for i = 1:N
    isNeighbour(i, :) = FindNeighbours(r, i, Rf, L);
end

% Exercise 8.8
%h = [1:25];
h = [25];
for H = 1:length(h)
    
oldThetas = {}; oldThetas{end + 1} = {theta};
for m = 1:S
    theta = UpdateOrientationWithDelay(theta, isNeighbour, eta, deltaT, h(H), oldThetas);
    oldThetas{end + 1} = {theta};
    velocities = UpdateVelocities(v, theta);
    deltaR = velocities.*deltaT;
    r = UpdatePositions(r, deltaR, L);
    globAlignCoeff(m) = CalculateGlobalAlignmentCoefficient(velocities, v);
    globClustCoeff(m) = CalculateGlobalClusteringCoefficent(r, Rf, L);
    
    isNeighbour = zeros(N) ~= 0;
    for i = 1:N
        isNeighbour(i, :) = FindNeighbours(r, i, Rf, L);
    end
end

subplot(1, 3, 1)
PlotVoronoiDiagram(initialR, L)
axis equal
axis([-3/2*L 3/2*L -3/2*L 3/2*L])
title('Initial distribution.')

subplot(1, 3, 2)
PlotVoronoiDiagram(r, L)
axis equal
axis([-3/2*L 3/2*L -3/2*L 3/2*L])
title('Distribution after 5*10^3 timesteps.')

subplot(1, 3, 3)
hold on
plot(globAlignCoeff)
plot(globClustCoeff)
legend('globAlignCoeff', 'globClustCoeff')
hold off

end

toc
