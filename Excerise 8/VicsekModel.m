%% Standard Vicsek model.
clear all
tic
% Setup.
N = 1000; % Number of particles. 
L = 100; % Side length. 
Rf = 10; % Interaction radius. 
eta = 0.1; % Noise level.
deltaT = 1; % Time step.
v = 1; % Speed.
rho = N/L^2; % Particle density.
S = 10^4; % Time steps.

% Exercise 8.1a)
%r = InitializePositions(N, L); % Random initialization. 
%save('initialR', 'r');
%save('initialR1000', 'r');
%load('initialR', 'r');
load('initialR1000', 'r');
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
%save('initialTheta1000', 'theta');
%load('initialTheta', 'theta');
load('initialTheta1000', 'theta');
globAlignCoeff = zeros(1, S);
globClustCoeff = zeros(1, S);

velocities = UpdateVelocities(v, theta);
deltaR = velocities.*deltaT;
r = UpdatePositions(r, deltaR, L);
theta = UpdateOrientation(theta, isNeighbour, eta, deltaT);
globAlignCoeff(1) = CalculateGlobalAlignmentCoefficient(velocities, v);
globClustCoeff(1) = CalculateGlobalClusteringCoefficent(r, Rf, L);

% Exercise 8.3
%PlotVoronoiDiagram(r, L)

%% Exercise 8.4a)
for m = 2:S
    isNeighbour = zeros(N) ~= 0;
    for i = 1:N
        isNeighbour(i, :) = FindNeighbours(r, i, Rf, L);
    end
    
    velocities = UpdateVelocities(v, theta);
    deltaR = velocities.*deltaT;
    r = UpdatePositions(r, deltaR, L);
    theta = UpdateOrientation(theta, isNeighbour, eta, deltaT);
    globAlignCoeff(m) = CalculateGlobalAlignmentCoefficient(velocities, v);
    globClustCoeff(m) = CalculateGlobalClusteringCoefficent(r, Rf, L);
    
    if (m == 10)
        r10 = r;
    end
    if (m == 100)
        r100 = r;
    end
    if (m == 500)
        r500 = r;
    end
    if (m == 1000)
        r1000 = r;
    end
    if (rem(m, 1000) == 0)
        disp(m)
    end
end

subplot(2, 5, 1)
PlotVoronoiDiagram(initialR, L)
axis equal
axis([-3/2*L 3/2*L -3/2*L 3/2*L])
title('Initial distribution.')

subplot(2, 5, 2)
PlotVoronoiDiagram(r10, L)
axis equal
axis([-3/2*L 3/2*L -3/2*L 3/2*L])
title('Distribution after 10 timesteps.')

subplot(2, 5, 3)
PlotVoronoiDiagram(r100, L)
axis equal
axis([-3/2*L 3/2*L -3/2*L 3/2*L])
title('Distribution after 100 timesteps.')

subplot(2, 5, 6)
PlotVoronoiDiagram(r500, L)
axis equal
axis([-3/2*L 3/2*L -3/2*L 3/2*L])
title('Distribution after 500 timesteps.')

subplot(2, 5, 7)
PlotVoronoiDiagram(r1000, L)
axis equal
axis([-3/2*L 3/2*L -3/2*L 3/2*L])
title('Distribution after 10^3 timesteps.')

subplot(2, 5, 8)
PlotVoronoiDiagram(r, L)
axis equal
axis([-3/2*L 3/2*L -3/2*L 3/2*L])
title('Distribution after 10^4 timesteps.')

subplot(2, 5, [4, 10])
hold on
plot(globAlignCoeff)
plot(globClustCoeff)
legend('ψ', 'c')
hold off

toc
