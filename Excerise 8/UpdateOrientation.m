function updatedTheta = UpdateOrientation(theta, isNeighbour, eta, deltaT)
    N = size(theta, 1);
    r = rand - 0.5;
    updatedTheta = zeros(N, 1);
    for i = 1:N
        neighbours = find(isNeighbour(i,:));
        averageTheta = atan2(mean(sin(theta(neighbours))), mean(cos(theta(neighbours))));
        updatedTheta(i) = averageTheta + eta*r*deltaT;
    end
end