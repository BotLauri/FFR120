%% Standard Vicsec model.
clear all
% Setup.
N = 10; % Number of particles. 
L = 10; % Side length. 
Rf = 3; % Interaction radius. 

% Exercise 8.1a)
particles = InitializePositions(N, L);
%deltaR = [1 1; 1 1; 1 1; 1 1; 1 1;];
%particles = UpdatePositions(particles, deltaR, L); 

% Exercise 8.1b)
index = 3; 
isNeighbour = FindNeighbours(particles, index, Rf, L);

% Exercise 8.2
globAlignCoeff = CalculationGlobalAlignmentCoefficient(particles, 
