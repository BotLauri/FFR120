%% Standard Vicsec model.
clear all
tic
% Setup.
N = 100; % Number of particles. 
L = 100; % Side length. 
Rf = 1; % Interaction radius. 
eta = 0.01; % Noise level.
deltaT = 1; % Time step.
v = 1; % Speed.
rho = N/L^2; % Particle density.
S = 10^3; % Time steps.

% Exercise 8.1a)
%r = InitializePositions(N, L); % Random initialization. 
%save('initialR', 'r');
load('initialR', 'r');
initialR = r;
%r = UpdatePositions(r, deltaR, L); 

% Exercise 8.1b)
isNeighbour = zeros(N) ~= 0;
for i = 1:N
    isNeighbour(i, :) = FindNeighbours(r, i, Rf, L);
end

% Exercise 8.2
%theta = InitializeOrientations(N); % Random initialization. 
%save('initialTheta', 'theta');
load('initialTheta', 'theta');
globAlignCoeff = zeros(1, S);
globClustCoeff = zeros(1, S);

theta = UpdateOrientation(theta, isNeighbour, eta, deltaT);
velocities = UpdateVelocities(v, theta);
deltaR = velocities.*deltaT;
r = UpdatePositions(r, deltaR, L);
globAlignCoeff(1) = CalculateGlobalAlignmentCoefficient(velocities, v);
globClustCoeff(1) = CalculateGlobalClusteringCoefficent(r, Rf, L);

% Exercise 8.3
%PlotVoronoiDiagram(r, L)

%% Exercise 8.4a)
for m = 2:S
    theta = UpdateOrientation(theta, isNeighbour, eta, deltaT);
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
title('Distribution after 10^4 timesteps.')

subplot(1, 3, 3)
hold on
plot(globAlignCoeff)
plot(globClustCoeff)
legend('globAlignCoeff', 'globClustCoeff')
hold off

toc
