function isNeighbour = FindNeighbours(particles, index, Rf, L)
    isNeighbour = zeros(1, size(particles, 1)) ~= 0;
    currentParticle = particles(index, :);
    % If the absolute distance is over L/2 then we have wrapped around. 
    for i = 1:size(particles, 1)
        xDistance = currentParticle(1) - particles(i, 1);
        if (abs(xDistance) > L/2)
            xDistance = L - xDistance;
        end
        yDistance = currentParticle(2) - particles(i, 2);
        if (abs(yDistance) > L/2)
            yDistance = L - yDistance;
        end
        distance = sqrt(xDistance^2 + yDistance^2);
        if (distance < Rf)
            isNeighbour(i) = true;
        end
    end
end

