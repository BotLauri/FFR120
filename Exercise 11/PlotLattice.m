function PlotLattice(lattice)
    positions = numel(lattice);
    cells = zeros(sqrt(positions), sqrt(positions));
    map = [1 1 1
           0 0 1
           1 0.647 0 % Orange
           0 1 0];
    for i = 1:positions
        mostFrequentAgent = mode(lattice{i});
        if (isnan(mostFrequentAgent))
            cells(i) = 0;
        elseif (mostFrequentAgent == 1)
            cells(i) = 1;
        elseif (mostFrequentAgent == 2)
            cells(i) = 2;
        elseif (mostFrequentAgent == 3)
            cells(i) = 3;
        end
    end
    colormap(map)
    imagesc(cells)
end

