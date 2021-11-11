function globalAlignmentCoefficient = CalculateGlobalAlignmentCoefficient(velocities, v)
    N = size(velocities, 1);
    term = zeros(N, 1);
    for j = 1:N
        term(j) = sqrt((velocities(j, 1)/v)^2 + (velocities(j, 2)/v)^2);
    end
    globalAlignmentCoefficient = 1/N * sum(term);
end

