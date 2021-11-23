%% SIR model with temporary immunity.
clear all
tic
L = 100; % Length of the lattice.
N = 1000; % Number of agents. 
initialInfectionRate = 0.01;
d = 0.8; % Probability of moving. 

% Parameters for 11.4a.
gammas = 0.01; % Probability of recovery.
betas = 0.2; % Probability of infection. 
alphas = 0.005; % Probability of losing immunity when infected.
params = [betas, gammas, alphas];
trials = 1;

% Parameters for 11.4b.
betas = [0.1, 0.2]; 
gammas = [0.01, 0.02]; 
alphas = [0.001:0.005:0.05]; 
[betaList, gammaList, alphaList] = meshgrid(betas, gammas, alphas);
params = [betaList(:), gammaList(:), alphaList(:)];
trials = 5;

for m = 1:trials
disp(m)

data = zeros(length(params), 7);
for T = 1:size(params, 1)
beta = params(T, 1); gamma = params(T, 2); alpha = params(T, 3);
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
    lattice = ImmunityLoss(lattice, alpha);
    SIR(t + 1, :) = GetData(lattice);
    if (rem(t, 10) == 0)
        %Plots(lattice, SIR, beta, gamma, t);
        %pause(0.01)
    end
    t = t + 1;
    if (rem(t, 1000) == 0)
        %disp(t)
    end 
end

SIR(t + 1, :) = GetData(lattice);
%Plots(lattice, SIR, beta, gamma, t);
%title('Example simulation with \alpha = 0.005, \beta = 0.2, \gamma = 0.01')

data(T, :) = [beta, gamma, alpha, SIR(end, 1), SIR(end, 2), SIR(end, 3), t];
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

%save('dataAlpha', 'avgData');

toc

%% R_inf as a function of mu. 
load('dataAlpha', 'avgData');
index11 = (avgData(:, 1) == 0.1).*(avgData(:, 2) == 0.01);
index12 = (avgData(:, 1) == 0.1).*(avgData(:, 2) == 0.02);
index21 = (avgData(:, 1) == 0.2).*(avgData(:, 2) == 0.01);
index22 = (avgData(:, 1) == 0.2).*(avgData(:, 2) == 0.02);
hold on 
scatter(avgData(find(index11), 3), avgData(find(index11), 6))
scatter(avgData(find(index12), 3), avgData(find(index12), 6))
scatter(avgData(find(index21), 3), avgData(find(index21), 6))
scatter(avgData(find(index22), 3), avgData(find(index22), 6))
title('R_\infty as a function of \mu averaged over 5 iterations for different \beta and \gamma.')
legend('\beta = 0.1, \gamma = 0.01', '\beta = 0.1, \gamma = 0.02', '\beta = 0.2, \gamma = 0.01', ... 
        '\beta = 0.2, \gamma = 0.02', 'Location', 'northeast')
xlabel('\mu')
ylabel('D_\infty')
hold off
