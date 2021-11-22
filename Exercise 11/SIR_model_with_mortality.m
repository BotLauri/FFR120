%% SIR model with mortality.
clear all
tic
L = 100; % Length of the lattice.
N = 1000; % Number of agents. 
initialInfectionRate = 0.01;
d = 0.8; % Probability of moving. 
%gammas = [0.01:0.02:0.1]; % Probability of recovery.
%betas = [0.01:0.1:1]; % Probability of infection. 
%mu = [0.001:0.02:0.02]; % Probability of death when infected. 
%[betaList, gammaList, muList] = meshgrid(betas, gammas, mu);
%params = [betaList(:), gammaList(:), muList(:)];
trials = 1;

for m = 1:trials

%data = zeros(length(gammas)*length(betas), 6);
%for T = 1:(length(gammas)*length(betas))
for T = 1:1
%beta = params(T, 1); gamma = params(T, 2);
beta = 0.2; gamma = 0.01; mu = 0.005;
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
    lattice = Death(lattice, gamma);
    SIR(t + 1, :) = GetData(lattice, mu);
    if (rem(t, 10) == 0)
        Plots(lattice, SIR, beta, gamma, t);
        pause(0.01)
    end
    t = t + 1;
end

SIR(t + 1, :) = GetData(lattice, mu);
Plots(lattice, SIR, beta, gamma, t);

data(T, :) = [beta, gamma, SIR(end, 1), SIR(end, 2), SIR(end, 3), t];
if (rem(T, 10) == 0)
    disp(T)
end

end % Loop through parameters. 

disp(trials)
if (trials == 1)
    avgData = data;
else
    avgData = (data(:, 3:7) + data(:, 3:7)) / 2;
end

end % Loop through trials. 

%save('data', 'data');

toc

%% R_inf as a function of beta. 
load('dataGamma_01', 'gamma01');
load('dataGamma_02', 'gamma02');
hold on 
scatter(data(:, 1), data(:, 5))
title('R_\infty as a function of \beta.')
xlabel('\beta')
ylabel('R_\infty')
hold off

%% Plot of parameters which give R_inf > 500.
isRinf500 = data(:, 5) > 500;
isRinf500 = reshape(isRinf500, [length(gammas), length(betas)]);
imagesc([betas(1) betas(end)], [gammas(1) gammas(end)], isRinf500)
xlabel('\beta')
ylabel('\gamma')
title('R_\infty > 500 as a function of \beta and \gamma.')
