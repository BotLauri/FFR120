%% SIR model.
clear all
tic
L = 100; % Length of the lattice.
N = 1000; % Number of agents. 
initialInfectionRate = 0.01;
d = 0.8; % Probability of moving. 
gammas = [0.01:0.007:0.1]; % Probability of recovery.
betas = [0.01:0.03:1]; % Probability of infection. 
[betaList, gammaList] = meshgrid(betas, gammas);
params = [betaList(:), gammaList(:)];
trials = 1;

betaGammas = [1:3:100];
[betaList, betaGammaList] = meshgrid(betas, betas./betaGammas);
params = [betaList(:), betaGammaList(:)];

%gammas = 0.02;
%betas = [0.01:0.02:1];
%trials = 5;

for m = 1:trials
disp(m)

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

%save('dataVaryingBetaGamma', 'data');
%save('dataVaryingBetaBetagamma', 'data');
%save('dataGamma_01', 'avgData');
%save('dataGamma_02', 'avgData');

toc

%% R_inf as a function of beta. 
load('dataGamma_01', 'avgData');
gamma01 = avgData;
load('dataGamma_02', 'avgData');
gamma02 = avgData;
hold on 
scatter(gamma01(:, 1), gamma01(:, 5))
scatter(gamma02(:, 1), gamma02(:, 5))
title('R_\infty as a function of \beta averaged over 5 iterations.')
legend('\gamma = 0.01', '\gamma = 0.02', 'Location', 'northwest')
xlabel('\beta')
ylabel('R_\infty')
hold off

%% Plot of parameters beta and gamma.
load('dataVaryingBetaGamma', 'data');
Rinf = data(:, 5);
Rinf = reshape(Rinf, [length(gammas), length(betas)]);
imagesc([betas(1) betas(end)], [gammas(1) gammas(end)], Rinf)
colorbar
xlabel('\beta')
ylabel('\gamma')
title('R_\infty as a function of \beta and \gamma.')

%% Plot of parameters beta and beta/gamma.
load('dataVaryingBetaBetagamma', 'data');
Rinf = data(:, 5);
Rinf = reshape(Rinf, [length(gammas), length(betas)]);
imagesc([betas(1) betas(end)], [gammas(1) gammas(end)], Rinf)
colorbar
xlabel('\beta')
ylabel('\beta / \gamma')
title('R_\infty as a function of \beta and \beta / \gamma.')
