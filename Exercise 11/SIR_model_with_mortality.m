%% SIR model with mortality.
clear all
tic
L = 100; % Length of the lattice.
N = 1000; % Number of agents. 
initialInfectionRate = 0.01;
d = 0.8; % Probability of moving. 

% Parameters for 11.3.
betas = [0.2, 0.5, 1]; 
gammas = [0.001, 0.01, 0.02]; 
mu = [0.0001:0.0007:0.02]; % Probability of death when infected.
[betaList, gammaList, muList] = meshgrid(betas, gammas, mu);
params = [betaList(:), gammaList(:), muList(:)];
trials = 3;

for m = 1:trials
disp(m)

data = zeros(length(params), 8);
for T = 1:length(params)
beta = params(T, 1); gamma = params(T, 2); mu = params(T, 3);
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
while (~CheckForInfected(lattice) && (t < 5000))
    lattice = Diffusion(lattice, d);
    lattice = Infection(lattice, beta);
    lattice = Recovery(lattice, gamma);
    lattice = Death(lattice, mu);
    SIR(t + 1, :) = GetData(lattice, mu);
    t = t + 1;
end

SIR(t + 1, :) = GetData(lattice, mu);

data(T, :) = [beta, gamma, mu, SIR(end, 1), SIR(end, 2), SIR(end, 3), SIR(end, 4), t];
if (rem(T, 10) == 0)
    disp(T)
end

end % Loop through parameters. 

if (m == 1)
    avgData = data;
else
    avgData = (avgData(:, :) + data(:, :)) / 2;
end

end % Loop through trials. 

save('dataMu', 'avgData');

toc

%% D_inf as a function of mu. 
load('dataMu', 'avgData');
index11 = (avgData(:, 1) == 0.2).*(avgData(:, 2) == 0.001);
index12 = (avgData(:, 1) == 0.2).*(avgData(:, 2) == 0.01);
index13 = (avgData(:, 1) == 0.2).*(avgData(:, 2) == 0.02);
index21 = (avgData(:, 1) == 0.5).*(avgData(:, 2) == 0.001);
index22 = (avgData(:, 1) == 0.5).*(avgData(:, 2) == 0.01);
index23 = (avgData(:, 1) == 0.5).*(avgData(:, 2) == 0.02);
index31 = (avgData(:, 1) == 1).*(avgData(:, 2) == 0.001);
index32 = (avgData(:, 1) == 1).*(avgData(:, 2) == 0.01);
index33 = (avgData(:, 1) == 1).*(avgData(:, 2) == 0.02);
hold on 
scatter(avgData(find(index11), 3), avgData(find(index11), 7))
scatter(avgData(find(index12), 3), avgData(find(index12), 7))
scatter(avgData(find(index13), 3), avgData(find(index13), 7))
scatter(avgData(find(index21), 3), avgData(find(index21), 7))
scatter(avgData(find(index22), 3), avgData(find(index22), 7))
scatter(avgData(find(index23), 3), avgData(find(index23), 7))
scatter(avgData(find(index31), 3), avgData(find(index31), 7))
scatter(avgData(find(index32), 3), avgData(find(index32), 7))
scatter(avgData(find(index33), 3), avgData(find(index33), 7))
title('D_\infty as a function of \mu for different \beta and \gamma averaged over 3 iterations.')
legend('\beta = 0.2, \gamma = 0.001', '\beta = 0.2, \gamma = 0.01', '\beta = 0.2, \gamma = 0.02', ... 
        '\beta = 0.5, \gamma = 0.001', '\beta = 0.5, \gamma = 0.01', '\beta = 0.5, \gamma = 0.02', ... 
         '\beta = 1, \gamma = 0.001',   '\beta = 1, \gamma = 0.01',   '\beta = 1, \gamma = 0.02', 'Location', 'northwest')
xlabel('\mu')
ylabel('D_\infty')
hold off
