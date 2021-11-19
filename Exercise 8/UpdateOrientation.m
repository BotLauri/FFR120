function updatedTheta = UpdateOrientation(theta, isNeighbour, eta, deltaT)
    N = size(theta, 1);
    updatedTheta = zeros(N, 1);
    for i = 1:N
        r = rand - 0.5;
        neighbours = find(isNeighbour(i,:));
        if (length(neighbours) == 1)
            averageTheta = theta(i);
        else
            averageTheta = atan2(mean(sin(theta(neighbours))), mean(cos(theta(neighbours))));
        end
        updatedTheta(i) = averageTheta + eta*r*deltaT;
    end
end