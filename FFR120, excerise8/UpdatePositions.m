function particles = UpdatePositions(deltaR, particles, L)
    for i = 1:size(particles, 1)
        particles(i, :) = particles(i, :) + deltaR(i, :);
    end
    for i = 1:size(particles, 1)
        for j = 1:size(particles, 2)
            if (particles(i, j) > L/2)
                particles(i, j) = particles(i, j) - L;
            elseif (particles(i ,j) < -L/2)
                particles(i, j) = particles(i, j) + L;
            end
        end
    end
end

