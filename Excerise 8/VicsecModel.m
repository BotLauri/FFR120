%% Standard Vicsec model.
clear all
tic
% Setup.
N = 10; % Number of particles. 
L = 10; % Side length. 
Rf = 3; % Interaction radius. 
eta = 0.1; % Noise level.
deltaT = 0.1; % Time step.
v = 1; % Speed.
rho = N/L^2; % Particle density.
n = 10; % Time steps.

% Exercise 8.1a)
r = InitializePositions(N, L); % Random initialization. 
%r = UpdatePositions(r, deltaR, L); 

% Exercise 8.1b)
isNeighbour = zeros(N) ~= 0;
for i = 1:N
    isNeighbour(i, :) = FindNeighbours(r, i, Rf, L);
end

% Exercise 8.2
theta = InitializeOrientations(N);
for m = 1:n
    theta = UpdateOrientation(theta, isNeighbour, eta, deltaT);
    velocities = UpdateVelocities(v, theta);
    deltaR = velocities.*deltaT;
    r = UpdatePositions(r, deltaR, L);
    globAlignCoeff = CalculateGlobalAlignmentCoefficient(velocities, v);
    %pause(1)
end

toc
