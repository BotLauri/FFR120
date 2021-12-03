%% Prisoner's dilemma on a lattice.
clear all
tic
L = 30; N = 7; mu = 0; timesteps = 20; R = 0.9;
numberOfInitialDefectors = 3*3;

% Initialization. 
if (numberOfInitialDefectors == 1)
    lattice = N*ones(L);
    lattice(ceil(L/2), ceil(L/2)) = 0;
    originalLattice = lattice;
end
if (numberOfInitialDefectors == 2)
    lattice = N*ones(L);
    lattice(ceil(L/3), ceil(2*L/3)) = 0;
    lattice(ceil(2*L/3), ceil(L/3)) = 0;
    originalLattice = lattice;
end
if (numberOfInitialDefectors == 3)
    lattice = N*ones(L);
    lattice(ceil(L/2), ceil(L/2)) = 0;
    lattice(ceil(L/4), ceil(3*L/4)) = 0;
    lattice(ceil(3*L/4), ceil(L/4)) = 0;
    originalLattice = lattice;
end
if (numberOfInitialDefectors == 4)
    lattice = N*ones(L);
    lattice(ceil(L/5), ceil(4*L/5)) = 0;
    lattice(ceil(2*L/5), ceil(3*L/5)) = 0;
    lattice(ceil(3*L/5), ceil(2*L/5)) = 0;
    lattice(ceil(4*L/5), ceil(L/5)) = 0;
    originalLattice = lattice;
end
if (numberOfInitialDefectors == 6)
    lattice = N*ones(L);
    lattice(ceil(L/5), ceil(4*L/5)) = 0;
    lattice(ceil(2*L/5), ceil(3*L/5)) = 0;
    lattice(ceil(3*L/5), ceil(2*L/5)) = 0;
    lattice(ceil(4*L/5), ceil(L/5)) = 0;
    lattice(ceil(L/5), ceil(L/5)) = 0;
    lattice(ceil(4*L/5), ceil(4*L/5)) = 0;
    originalLattice = lattice;
end
if (numberOfInitialDefectors == L*L-1)
    lattice = zeros(L);
    lattice(ceil(L/2), ceil(L/2)) = N;
    originalLattice = lattice;
end
if (numberOfInitialDefectors == 3*3)
    lattice = zeros(L);
    lattice(ceil(L/2)-1, ceil(L/2)-1) = N;
    lattice(ceil(L/2)-1, ceil(L/2)) = N;
    lattice(ceil(L/2)-1, ceil(L/2)+1) = N;
    lattice(ceil(L/2), ceil(L/2)-1) = N;
    lattice(ceil(L/2), ceil(L/2)) = N;
    lattice(ceil(L/2), ceil(L/2)+1) = N;
    lattice(ceil(L/2)+1, ceil(L/2)-1) = N;
    lattice(ceil(L/2)+1, ceil(L/2)) = N;
    lattice(ceil(L/2)+1, ceil(L/2)+1) = N;
    originalLattice = lattice;
end

for t = 1:timesteps
    
    % Competition.
    competitionLattice = zeros(L);
    for i = 1:L
        for j = 1:L
            if (j - 1 == 0)
                west = PrisonersDilemma(lattice(i, j), lattice(i, L), N, R);
            else
                west = PrisonersDilemma(lattice(i, j), lattice(i, j - 1), N, R);
            end
            if (j + 1 == L + 1)
                east = PrisonersDilemma(lattice(i, j), lattice(i, 1), N, R);
            else
                east = PrisonersDilemma(lattice(i, j), lattice(i, j + 1), N, R);
            end
            if (i + 1 == L + 1)
                south = PrisonersDilemma(lattice(i, j), lattice(1, j), N, R);
            else
                south = PrisonersDilemma(lattice(i, j), lattice(i + 1, j), N, R);
            end
            if (i - 1 == 0)
                north = PrisonersDilemma(lattice(i, j), lattice(L, j), N, R);
            else
                north = PrisonersDilemma(lattice(i, j), lattice(i - 1, j), N, R);
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
    
end % Time loop.

colormap(flipud(gray))
imagesc(lattice)

toc

%% Plots.
colormap(flipud(gray))
subplot(1, 2, 1)
imagesc(originalLattice)
colorbar
ylabel('t = 0')
title('Initial defectors: ', strcat(num2str(numberOfInitialDefectors, 2)))
pbaspect([1 1 1])

subplot(1, 2, 2)
imagesc(lattice)
colorbar
ylabel('t = 20')
title('Initial defectors: ', strcat(num2str(numberOfInitialDefectors, 2)))
pbaspect([1 1 1])

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
ylabel('t = 30')
title('R: ', strcat(num2str(R, 2)))
pbaspect([1 1 1])
