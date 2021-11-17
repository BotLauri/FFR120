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
T = 1; % Trials.

for t = 1:T

%r = InitializePositions(N, L); % Random initialization. 
%save('initialRDelay', 'r');
load('initialRDelay', 'r');
initialR = r;
%theta = InitializeOrientations(N); % Random initialization. 
load('initialTheta', 'theta');

% Exercise 8.8 / 8.9
H = [0:25]; % Change this depending on what excerise (8.8 / 8.9).
%H = [-15:-1]; % Remember to change UpdateOrientationWithDelay aswell!
globAlignCoeff = zeros(length(H), S);
globClustCoeff = zeros(length(H), S);

for h = 1:length(H)
%load('initialRDelay', 'r');    
%load('initialTheta', 'theta');
    
oldThetas = {}; oldThetas{end + 1} = {theta};
for m = 1:S
    isNeighbour = zeros(N) ~= 0;
    for i = 1:N
        isNeighbour(i, :) = FindNeighbours(r, i, Rf, L);
    end
    
    if (H(h) < 0)
        predTheta = UpdateOrientationWithDelay(theta, isNeighbour, eta, deltaT, H(h), oldThetas);
        velocities = UpdateVelocities(v, predTheta);
        deltaR = velocities.*deltaT;
        r = UpdatePositions(r, deltaR, L);
        theta = UpdateOrientation(theta, isNeighbour, eta, deltaT);
    else
        velocities = UpdateVelocities(v, theta);
        deltaR = velocities.*deltaT;
        r = UpdatePositions(r, deltaR, L);
        theta = UpdateOrientationWithDelay(theta, isNeighbour, eta, deltaT, H(h), oldThetas);
    end
    oldThetas{end + 1} = {theta};
    globAlignCoeff(h, m) = CalculateGlobalAlignmentCoefficient(velocities, v);
    globClustCoeff(h, m) = CalculateGlobalClusteringCoefficent(r, Rf, L);
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

if (H(h) == -1)
    rhm1 = r;
end

if (H(h) == -5)
    rhm5 = r;
end

if (H(h) == -15)
    rhm15 = r;
end

disp(H(h))

end % H loop.

if (t == 1)
    oldGlobAlignCoeff = globAlignCoeff;
    oldGlobClustCoeff = globClustCoeff;
else
    meanGlobAlignCoeff = (oldGlobAlignCoeff + globAlignCoeff)/2;
    meanGlobClustCoeff = (oldGlobClustCoeff + globClustCoeff)/2;
    oldGlobAlignCoeff = globAlignCoeff;
    oldGlobClustCoeff = globClustCoeff;
end

t

end % T loop.

toc

%% Plots positive h.

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

%% Plots negative h.

subplot(1, 4, 1)
PlotVoronoiDiagram(initialR, L)
axis equal
axis([-3/2*L 3/2*L -3/2*L 3/2*L])
title('Initial distribution.')

subplot(1, 4, 2)
PlotVoronoiDiagram(rhm1, L)
axis equal
axis([-3/2*L 3/2*L -3/2*L 3/2*L])
title('Distribution after 10^4 timesteps for h = -1.')

subplot(1, 4, 3)
PlotVoronoiDiagram(rhm5, L)
axis equal
axis([-3/2*L 3/2*L -3/2*L 3/2*L])
title('Distribution after 10^4 timesteps for h = -5.')

subplot(1, 4, 4)
PlotVoronoiDiagram(rhm15, L)
axis equal
axis([-3/2*L 3/2*L -3/2*L 3/2*L])
title('Distribution after 10^4 timesteps for h = -15.')

%save('globalAlignmentCoefficientMinus', 'globAlignCoeff');
%save('globalClusteringCoefficientMinus', 'globClustCoeff');

%% Plot statistics. 

globAlignCoeffEnd = meanGlobAlignCoeff(:, 5001:10000);
globClustCoeffEnd = meanGlobClustCoeff(:, 5001:10000);
for i = 1:size(globClustCoeffEnd, 1)
    meanAlign(i) = mean(globAlignCoeffEnd(i, :));
    meanClust(i) = mean(globClustCoeffEnd(i, :));
    errorAlign(i) = std(globAlignCoeffEnd(i, :));
    errorClust(i) = std(globClustCoeffEnd(i, :));
end
hold on
errorbar(H, meanAlign, errorAlign);
errorbar(H, meanClust, errorClust);
legend('ψ', 'c')
hold off

%save('meanGlobalAlignmentCoefficient', 'meanGlobAlignCoeff');
%save('meanGlobalClusteringCoefficient', 'meanGlobClustCoeff');

%% Other plots.

hold on
plot(globAlignCoeff(3, :))
plot(globClustCoeff(3, :))
legend('ψ', 'c')
hold off
