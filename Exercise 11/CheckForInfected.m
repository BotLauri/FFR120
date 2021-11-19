function endSimulation = CheckForInfected(lattice)
    endSimulation = true;
    L = length(lattice);
    for i = 1:L
        for j = 1:L
            if (~isempty(find(lattice{i,j} == 2, 1)))
                endSimulation = false;
                return
            end
        end
    end
end