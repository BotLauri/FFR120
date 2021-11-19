%% Vicsec model with delay.
clear all
tic
% Setup.
N = 100; % Number of particles. 
L = 100; % Side length. 
Rf = 1; % Interaction radius.
eta = 0.1; % Noise level.
deltaT = 1; % Time step.
v = 1; % Speed.
rho = N/L^2; % Particle density.
S = 10^4; % Time steps.
K = [1, 4, 8]; % Nearest neighbours. 
vision = 50; % Vision in degrees.

load('initialR', 'r');
initialR = r;
load('initialTheta', 'theta');

globAlignCoeff = zeros(length(K), S);
globClustCoeff = zeros(length(K), S);
for k = 1:length(K)

for m = 1:S
    
    velocities = UpdateVelocities(v, theta);
    deltaR = velocities.*deltaT;
    r = UpdatePositions(r, deltaR, L);
    theta = UpdateOrientationNearestNeighbours(theta, r, K(k), eta, deltaT, vision);
    globAlignCoeff(k, m) = CalculateGlobalAlignmentCoefficient(velocities, v);
    globClustCoeff(k, m) = CalculateGlobalClusteringCoefficent(r, Rf, L);
    
    if (rem(m, 1000) == 0)
        disp(m)
    end
end

if (K(k) == 1)
    rk1 = r;
end

if (K(k) == 4)
    rk4 = r;
end

if (K(k) == 8)
    rk8 = r;
end

end

toc

%% Other plots.

subplot(2, 2, 1)
PlotVoronoiDiagram(initialR, L)
axis equal
axis([-3/2*L 3/2*L -3/2*L 3/2*L])
title('Initial distribution.')

subplot(2, 2, 2)
PlotVoronoiDiagram(rk1, L)
axis equal
axis([-3/2*L 3/2*L -3/2*L 3/2*L])
title('Distribution after 10^4 timesteps for k = 1 with vision = 50°.')

subplot(2, 2, 3)
PlotVoronoiDiagram(rk4, L)
axis equal
axis([-3/2*L 3/2*L -3/2*L 3/2*L])
title('Distribution after 10^4 timesteps for k = 4 with vision = 50°.')

subplot(2, 2, 4)
PlotVoronoiDiagram(rk8, L)
axis equal
axis([-3/2*L 3/2*L -3/2*L 3/2*L])
title('Distribution after 10^4 timesteps for k = 8 with vision = 50°.')

%save('globalAlignmentCoefficientNonMetric', 'globAlignCoeff');
%save('globalClusteringCoefficientNonMetric', 'globClustCoeff');
