function lattice = Death(lattice, mu)
    L = length(lattice);
     for i = 1:L
         for j = 1:L
             if (~isempty(find(lattice{i,j} == 2, 1)))
                infectedAgentIndex = find(lattice{i,j} == 2);
                for n = infectedAgentIndex
                    r = rand();
                    if (mu > r)
                        lattice{i,j}(n) = 4;
                    end
                end
            end
        end
    end
end

