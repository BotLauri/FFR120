function updatedTheta = UpdateOrientationWithDelay(theta, isNeighbour, eta, deltaT, h, oldThetas)
    n = size(oldThetas, 2);
    if (h >= n)
        updatedTheta = theta;
        return
    end
    N = size(theta, 1);
    r = rand - 0.5;
    updatedTheta = zeros(N, 1);
    for i = 1:N
        neighbours = find(isNeighbour(i,:));
        averageTheta = atan2(mean(sin(oldThetas{1, n-h}{1}(neighbours))), ...
                             mean(cos(oldThetas{1, n-h}{1}(neighbours))));
        updatedTheta(i) = averageTheta + eta*r*deltaT;
    end
end