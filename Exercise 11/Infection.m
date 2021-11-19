function lattice = Infection(lattice, beta)
     L = length(lattice);
     for i = 1:L
         for j = 1:L
             r = rand();
             if (beta > r && ~isempty(find(lattice{i,j} == 2, 1)))
                 numberOfAgents = length(lattice{i,j});
                 for n = 1:numberOfAgents
                     if (lattice{i,j}(n) == 1)
                         lattice{i,j}(n) = 2;
                     end
                 end
             end
         end
     end
end
