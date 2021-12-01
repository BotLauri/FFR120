%% Prisoner's dilemma on a lattice.
clear all
tic
L = 7;
N = 7; T = 0; R = 0.9; P = 1; S = 1.5; mu = 0;
timesteps = 20;

lattice = zeros(L);
lattice(ceil(L/2), ceil(L/2)) = N;
% Initialization.
%for i = 1:L
%    for j = 1:L
%        lattice(i, j) = randi(N+1) - 1;
%    end
%end

%imagesc(lattice)
%pause(1)

for t = 1:timesteps
    
    % Competition.
    updatedLattice = lattice;
    for i = 1:L
        for j = 1:L
            if (j - 1 == 0)
                west = PrisonerDilemma(lattice(i, j), lattice(i, L), N);
            else
                west = PrisonerDilemma(lattice(i, j), lattice(i, j - 1), N);
            end
            if (j + 1 == L + 1)
                east = PrisonerDilemma(lattice(i, j), lattice(i, 1), N);
            else
                east = PrisonerDilemma(lattice(i, j), lattice(i, j + 1), N);
            end
            if (i + 1 == L + 1)
                south = PrisonerDilemma(lattice(i, j), lattice(1, j), N);
            else
                south = PrisonerDilemma(lattice(i, j), lattice(i + 1, j), N);
            end
            if (i - 1 == 0)
                north = PrisonerDilemma(lattice(i, j), lattice(L, j), N);
            else
                north = PrisonerDilemma(lattice(i, j), lattice(i - 1, j), N);
            end
            self =  PrisonerDilemma(lattice(i, j), lattice(i, j), N);
            index = find(imregionalmin([west, east, south, north, self]));
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
        %[west, east, south, north, self]    
        end
    end
    lattice = updatedLattice; % Revision.
    
    % Mutation.
    for i = 1:L
        for j = 1:L
            if (rand() < mu)
                lattice(i, j) = randi(N+1) - 1;
            end
        end
    end
    
end % Time loop.

imagesc(lattice)

toc
