%% Prisoner's dilemma on a lattice.
% 35-45 sec per simulation with parfor. 

clear all
tic
L = 30; N = 7; mu = 0.01; timesteps = 500;

Rs = 0.01:0.03:0.99;
Ss = 1:0.06:3; 
%Rs = [0.1, 0.3];
%Ss = [1, 2];

dataStrategies = zeros(timesteps,  N + 4, length(Rs)*length(Ss));

for r = 1:length(Rs)
    
disp(r)
dataCurrentStrategies = zeros(timesteps,  N + 4, length(Ss));
    
parfor s = 1:length(Ss)
    
R = Rs(r);
S = Ss(s);
% Initialization.
lattice = zeros(L);
for i = 1:L
    for j = 1:L
        lattice(i, j) = randi(N+1) - 1;
    end
end

for t = 1:timesteps
    
    % Competition.
    competitionLattice = zeros(L);
    for i = 1:L
        for j = 1:L
            if (j - 1 == 0)
                west = PrisonersDilemma(lattice(i, j), lattice(i, L), N, R, S);
            else
                west = PrisonersDilemma(lattice(i, j), lattice(i, j - 1), N, R, S);
            end
            if (j + 1 == L + 1)
                east = PrisonersDilemma(lattice(i, j), lattice(i, 1), N, R, S);
            else
                east = PrisonersDilemma(lattice(i, j), lattice(i, j + 1), N, R, S);
            end
            if (i + 1 == L + 1)
                south = PrisonersDilemma(lattice(i, j), lattice(1, j), N, R, S);
            else
                south = PrisonersDilemma(lattice(i, j), lattice(i + 1, j), N, R, S);
            end
            if (i - 1 == 0)
                north = PrisonersDilemma(lattice(i, j), lattice(L, j), N, R, S);
            else
                north = PrisonersDilemma(lattice(i, j), lattice(i - 1, j), N, R, S);
            end
            competitionLattice(i, j) = west + east + south + north;
        end
    end
    
    % Revision.
    updatedLattice = lattice;
    for i = 1:L
        for j = 1:L
            if (j == 1)
                westAvg = competitionLattice(i, L);
            else
                westAvg = competitionLattice(i, j-1);
            end
            if (j == L)
                eastAvg = competitionLattice(i, 1);
            else
                eastAvg = competitionLattice(i, j+1);
            end
            if (i == L)
                southAvg = competitionLattice(1, j);
            else
                southAvg = competitionLattice(i+1, j);
            end
            if (i == 1)
                northAvg = competitionLattice(L, j);
            else
                northAvg = competitionLattice(i-1, j);
            end
            selfAvg = competitionLattice(i, j);
            [~, index] = min([westAvg, eastAvg, southAvg, northAvg, selfAvg], [], 'all', 'linear');
            if (length(index) ~= 1)
                I = randsample([index], 1);
            else
                I = index;
            end
            if (I == 1)
                if (j - 1 == 0)
                    updatedLattice(i, j) = lattice(i, L);
                else
                    updatedLattice(i, j) = lattice(i, j - 1);
                end
            elseif (I == 2)
                if (j + 1 == L + 1)
                    updatedLattice(i, j) = lattice(i, 1);
                else
                    updatedLattice(i, j) = lattice(i, j + 1);
                end
            elseif (I == 3)
                if (i + 1 == L + 1)
                    updatedLattice(i, j) = lattice(1, j);
                else
                    updatedLattice(i, j) = lattice(i + 1, j);
                end
            elseif (I == 4)
                if (i - 1 == 0)
                    updatedLattice(i, j) = lattice(L, j);
                else
                    updatedLattice(i, j) = lattice(i - 1, j);
                end
            else
                updatedLattice(i, j) = lattice(i, j);
            end
        end
    end
    lattice = updatedLattice;
    
    % Mutation.
    for i = 1:L
        for j = 1:L
            if (rand() < mu)
                lattice(i, j) = randi(N+1) - 1;
            end
        end
    end
    
    currentStrategies = zeros(1, N+4);
    currentStrategies(1) = R;
    currentStrategies(2) = S;
    currentStrategies(3) = t;
    for i = 1:N+1
        currentStrategies(i+3) = sum(lattice(:) == i-1);
    end
    dataCurrentStrategies(t, :, s) = currentStrategies;
    
    if rem(t, 100) == 0
        disp(t)
    end
    
end % Time loop.
end % S loop. 

dataStrategies(:, :, (r - 1)*length(Ss) + 1:(r - 1)*length(Ss) + length(Ss)) = dataCurrentStrategies;

end % R loop.

%save('dataStrategies', 'dataStrategies');

toc

%% Plot for each of the variances. 
dataStrategies = struct2cell(load('dataStrategies', 'dataStrategies'));
data = dataStrategies{1};
data = data(101:500, :, :);
variances = var(data);
Rs = 0.01:0.03:0.99;
Ss = 1:0.06:3; 

hold on
for n = 1:N+1
    subplot(2, 4, n)
    varI = reshape(variances(:, n+3, :), [length(Ss), length(Rs)]);
    imagesc([0.01, 0.99], [1, 3], varI);
    set(gca, 'YDir', 'normal')
    colormap('jet')
    pbaspect([1 1 1])
    colorbar
    xlabel('R')
    ylabel('S')
    title(strcat("\sigma_", int2str(n-1), "^2"))
end 
hold off

%% Plot for sum_n > 10000. 
dataStrategies = struct2cell(load('dataStrategies', 'dataStrategies'));
data = dataStrategies{1};
data = data(101:500, :, :);
variances = var(data);

varSum = zeros(1, length(Ss)*length(Rs));
for t = 1:length(Ss)*length(Rs)
    varSum(t) = sum(variances(:, 4:11, t)) > 10000;
end
varSums = reshape(varSum, [length(Ss), length(Rs)]);
imagesc([0.01, 0.99], [1, 3], varSums);
set(gca, 'YDir', 'normal')
colormap(flipud(gray))
pbaspect([1 1 1])
colorbar
xlabel('R')
ylabel('S')
title("\Sigma\sigma_n^2 > 10000")
