%% SIR model.
clear all
tic
L = 100; % Length of the lattice.
N = 1000; % Number of agents. 
initialInfectionRate = 0.01;
d = 0.8; % Probability of moving. 
gammas = [0.01:0.01:0.1]; % Probability of recovery.
betas = [0.05:0.05:1]; % Probability of infection. 
[betaList, gammaList] = meshgrid(betas, gammas);
params = [betaList(:), gammaList(:)];

data = zeros(length(gammas)*length(betas), 6);
for T = 1:(length(gammas)*length(betas))
beta = params(T, 1); gamma = params(T, 2);
t = 0; % Number of timesteps. 
SIR = [];
    
% Initialization of the grid. 
lattice = cell(L);
for n = 1:(N - N*initialInfectionRate)
    xCoord = randi(L);
    yCoord = randi(L);
    lattice{xCoord, yCoord} = [lattice{xCoord, yCoord}, 1];
    if (n > (N - 2*N*initialInfectionRate))
        lattice{xCoord, yCoord} = [lattice{xCoord, yCoord}, 2];
    end
end

% Main loop.
while (~CheckForInfected(lattice))
    lattice = Diffusion(lattice, d);
    lattice = Infection(lattice, beta);
    lattice = Recovery(lattice, gamma);
    SIR(t + 1, :) = GetData(lattice);
    if (rem(t, 10) == 0)
        %Plots(lattice, SIR, beta, gamma, t);
        %pause(0.01)
    end
    t = t + 1;
end

SIR(t + 1, :) = GetData(lattice);
%Plots(lattice, SIR, beta, gamma, t);

data(T, :) = [beta, gamma, SIR(end, 1), SIR(end, 2), SIR(end, 3), t];
save('data', 'data');

end
        
toc
