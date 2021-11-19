function PlotVoronoiDiagram(r, L)
R = [[r(:, 1) - L, r(:, 2) + L]; [r(:, 1), r(:, 2) + L]; [r(:, 1) + L, r(:, 2) + L];
     [r(:, 1) - L, r(:, 2)];     [r(:, 1), r(:, 2)];     [r(:, 1) + L, r(:, 2)];
     [r(:, 1) - L, r(:, 2) - L]; [r(:, 1), r(:, 2) - L]; [r(:, 1) + L, r(:, 2) - L]];
hold on
[vertices, cells] = voronoin(R);
ignoredIndexes = zeros(1, size(R, 1));
for i = 1:size(R, 1)
    if (max(abs(R(i, 1))) < 1/2*L)
        if (max(abs(R(i, 2))) < 1/2*L)
            ignoredIndexes(i) = 1;
        end
    end
end
voronoi(R(:, 1), R(:, 2));
for i = 1:length(cells)
    if (ignoredIndexes(i) == 1)
        V1 = vertices(cells{i}, 1); 
        V2 = vertices(cells{i}, 2);
        patch(V1, V2, 2)
    end
end
axis([-3/2*L 3/2*L -3/2*L 3/2*L])
line([L/2, L/2], [3*L/2, -3*L/2], 'Color', 'red', 'LineStyle', '--')
line([-L/2, -L/2], [3*L/2, -3*L/2], 'Color', 'red', 'LineStyle', '--')
line([3*L/2, -3*L/2], [L/2, L/2], 'Color', 'red', 'LineStyle', '--')
line([3*L/2, -3*L/2], [-L/2, -L/2], 'Color', 'red', 'LineStyle', '--')
hold off
end

