%% Prisoner's dilemma on a lattice.
clear all
tic
L = 30; N = 7; mu = 0.01; timesteps = 100; R = 0.84; S = 3;

% Initialization.
lattice = zeros(L);
for i = 1:L
    for j = 1:L
        if (rand() < 0.5)
            lattice(i, j) = 0;
        else
            lattice(i, j) = N;
        end
    end
end
originalLattice = lattice;

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
                if (rand() < 0.5)
                    lattice(i, j) = 0;
                else
                    lattice(i, j) = N;
                end
            end
        end
    end
    
end % Time loop.

colormap(flipud(gray))
imagesc(lattice)

toc

%% Plots with different R.
colormap(flipud(gray))
subplot(1, 2, 1)
imagesc(originalLattice)
colorbar
ylabel('t = 0')
title('R: ', strcat(num2str(R, 2)))
pbaspect([1 1 1])

subplot(1, 2, 2)
imagesc(lattice)
colorbar
ylabel('t = 100')
title('R: ', strcat(num2str(R, 2)))
pbaspect([1 1 1])

%% Plots with different S (and R = 0.84).
colormap(flipud(gray))
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
