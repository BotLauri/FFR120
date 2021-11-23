%% SIR model.
clear all
tic
L = 100; % Length of the lattice.
N = 1000; % Number of agents. 
initialInfectionRate = 0.01;
d = 0.8; % Probability of moving. 

% Plotting parameters. 
%gammas = 0.01; % Probability of recovery.
%betas = 0.2; % Probability of infection. 
%params = [betas, gammas];
%trials = 1;

% Parameters for 11.1b.
%gammas = 0.02;
%betas = [0.01:0.02:1];
%trials = 5;

% Parameters for 11.1d.
%gammas = [0.01:0.007:0.1];
%betas = [0.01:0.03:1]; 
%[betaList, gammaList] = meshgrid(betas, gammas);
%params = [betaList(:), gammaList(:)];
%trials = 1;

% Parameters for 11.2c.
betaGammas = [1:3:80];
betas = [0.01:0.04:1];
params = zeros(length(betaGammas)*length(betas), 2);
for i = 1:length(betaGammas)
    params((((i-1)*length(betas) + 1):(i-1)*length(betas) + length(betas)), 1) = betas;
    params((((i-1)*length(betas) + 1):(i-1)*length(betas) + length(betas)), 2) = betas(:) ./ betaGammas(i);
end
trials = 2;

for m = 1:trials
disp(m)

data = zeros(length(params), 6);
for T = 1:length(params)
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
while (~CheckForInfected(lattice) && (t < 5000))
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
%save('dataVaryingBetaBetagamma', 'avgData');
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

%% R_inf as a function of beta/gamma.
load('dataGamma_01', 'avgData');
gamma01 = avgData;
load('dataGamma_02', 'avgData');
gamma02 = avgData;
hold on 
scatter(gamma01(:, 1) / gamma01(1, 2), gamma01(:, 5))
scatter(gamma02(:, 1) / gamma02(1, 2), gamma02(:, 5))
title('R_\infty as a function of \beta / \gamma averaged over 5 iterations.')
legend('\gamma = 0.01', '\gamma = 0.02', 'Location', 'northwest')
xlabel('\beta / \gamma')
ylabel('R_\infty')
hold off

%% Phase diagram of R_inf as a function of beta and beta.
lenGammas = 13; lenBetas = 34;
load('dataVaryingBetaGamma', 'data');
betas = data(:, 1); gammas = data(:, 2); Rinf = data(:, 5);
Rinf = reshape(Rinf, [lenGammas, lenBetas]);
imagesc([betas(1) betas(end)], [gammas(1) gammas(end)], Rinf)
colorbar
xlabel('\beta')
ylabel('\gamma')
title('R_\infty as a function of \beta and \gamma.')
set(gca,'YDir','normal')

%% Phase diagram of R_inf as a function of beta and beta/gamma.
lenBetaGammas = 27; lenBetas = 25;
load('dataVaryingBetaBetagamma', 'avgData');
betas = avgData(:, 1);
betaGammas = avgData(:, 1) ./ avgData(:, 2);
Rinf = avgData(:, 5);
Rinf = reshape(Rinf, [lenBetas, lenBetaGammas]);
imagesc([betas(1) betas(end)], [betaGammas(1) betaGammas(end)], Rinf')
colorbar
xlabel('\beta')
ylabel('\beta / \gamma')
title('R_\infty as a function of \beta and \beta / \gamma averaged over 2 iterations.')
set(gca,'YDir','normal')
