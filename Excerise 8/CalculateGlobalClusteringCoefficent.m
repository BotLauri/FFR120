function globalClusteringCoefficient = CalculateGlobalClusteringCoefficent(r, Rf, L)
R = [[r(:, 1) - L, r(:, 2) + L]; [r(:, 1), r(:, 2) + L]; [r(:, 1) + L, r(:, 2) + L];
     [r(:, 1) - L, r(:, 2)];     [r(:, 1), r(:, 2)];     [r(:, 1) + L, r(:, 2)];
     [r(:, 1) - L, r(:, 2) - L]; [r(:, 1), r(:, 2) - L]; [r(:, 1) + L, r(:, 2) - L]];
[vertices, cells] = voronoin(R);
ignoredIndexes = zeros(1, size(R, 1));
for i = 1:size(R, 1)
    if (max(abs(R(i, 1))) < 1/2*L)
        if (max(abs(R(i, 2))) < 1/2*L)
            ignoredIndexes(i) = 1;
        end
    end
end
A = zeros(sum(ignoredIndexes), 1);
index = 1;
for i = find(ignoredIndexes)
    V1 = vertices(cells{i}, 1); 
    V2 = vertices(cells{i}, 2);
    A(index) = polyarea(V1, V2);
    index = index + 1;
end
globalClusteringCoefficient = sum(A < pi*Rf^2)/length(A);
end