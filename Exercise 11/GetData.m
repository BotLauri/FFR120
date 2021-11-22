function SIR = GetData(lattice, mu)
    if ~exist('mu', 'var')
    SIR = [0, 0, 0];
    L = length(lattice);
    for i = 1:L
        for j = 1:L
            numberOfAgents = size(lattice{i,j}, 2);
            for n = 1:numberOfAgents
                if (lattice{i,j}(n) == 1)
                    SIR(1) = SIR(1) + 1;
                elseif (lattice{i,j}(n) == 2)
                    SIR(2) = SIR(2) + 1;
                else
                    SIR(3) = SIR(3) + 1;
                end
            end
        end
    end
    else
    SIR = [0, 0, 0, 0];
    L = length(lattice);
    for i = 1:L
        for j = 1:L
            numberOfAgents = size(lattice{i,j}, 2);
            for n = 1:numberOfAgents
                if (lattice{i,j}(n) == 1)
                    SIR(1) = SIR(1) + 1;
                elseif (lattice{i,j}(n) == 2)
                    SIR(2) = SIR(2) + 1;
                elseif (lattice{i,j}(n) == 3)
                    SIR(3) = SIR(3) + 1;
                else
                    SIR(4) = SIR(4) + 1;
                end
            end
        end
    end
    end
end

