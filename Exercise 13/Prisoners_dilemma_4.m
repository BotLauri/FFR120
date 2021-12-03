%% Prisoner's dilemma on a lattice.
clear all
tic
L = 30; N = 7; mu = 0.01; timesteps = 100; R = 0.8; S = 1.5;

% Initialization.
lattice = zeros(L);
for i = 1:L
    for j = 1:L
        lattice(i, j) = randi(N+1) - 1;
    end
end

dataStrategies = zeros(timesteps, N+1);
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
    
    currentStrategies = zeros(1, N+1);
    for i = 1:N+1
        currentStrategies(i) = sum(lattice(:) == i-1);
    end
        
    dataStrategies(t, :) = currentStrategies;
    
end % Time loop.

colormap('jet')
imagesc(lattice)

toc

%% Plots with different R.
colormap('jet')
subplot(1, 2, 1)
imagesc(lattice)
colorbar
ylabel('t = 100')
title('R: ', strcat(num2str(R, 2)))
pbaspect([1 1 1])

colors = jet(N+1);
subplot(1, 2, 2)
hold on
for i = 1:N+1
    plot(1:1:timesteps, dataStrategies(:,i), 'color', colors(i,:))
end
legend('n = 0', 'n = 1', 'n = 2', 'n = 3', 'n = 4','n = 5', 'n = 6', 'n = 7')
xlabel('t')
ylabel('N_n')
title('Number of agents (N) with a particular strategy (n).')
hold off
pbaspect([1 1 1])

%% Plots with different S (and R = 0.84).
colormap('jet')
subplot(1, 2, 1)
imagesc(originalLattice)
colorbar
ylabel('t = 0')
title('R = 0.84, S: ', strcat(num2str(S, 3)))
pbaspect([1 1 1])

subplot(1, 2, 2)
imagesc(lattice)
colorbar
ylabel('t = 100')
title('R = 0.84, S: ', strcat(num2str(S, 3)))
pbaspect([1 1 1])
