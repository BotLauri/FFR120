function globalAlignmentCoefficient = CalculateGlobalAlignmentCoefficient(velocities, v)
    N = size(velocities, 1);
    globalAlignmentCoefficient = norm(sum(velocities) / v) / N;
end

