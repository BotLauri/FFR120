function lattice = Recovery(lattice, gamma)
    L = length(lattice);
     for i = 1:L
         for j = 1:L
             r = rand();
             if (gamma > r && ~isempty(find(lattice{i,j} == 2, 1)))
                 infectedAgentIndex = find(lattice{i,j} == 2);
                 for n = infectedAgentIndex
                     lattice{i,j}(n) = 3;
                 end
             end
         end
     end
end
