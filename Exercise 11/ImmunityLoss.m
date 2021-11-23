function lattice = ImmunityLoss(lattice, alpha)
    L = length(lattice);
     for i = 1:L
         for j = 1:L
             if (~isempty(find(lattice{i,j} == 3, 1)))
                infectedAgentIndex = find(lattice{i,j} == 3);
                for n = infectedAgentIndex
                    r = rand();
                    if (alpha > r)
                        lattice{i,j}(n) = 1;
                    end
                end
            end
        end
    end
end

