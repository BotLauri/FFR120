function updatedTheta = UpdateOrientationWithDelay(theta, isNeighbour, eta, deltaT, h, oldThetas)
    n = size(oldThetas, 2);
    if (h >= n)
        updatedTheta = theta;
        return
    end
    N = size(theta, 1);
    updatedTheta = zeros(N, 1);
    if (h <= 0)
        for i = 1:N
            r = rand - 0.5;
            neighbours = find(isNeighbour(i,:));
            if (length(neighbours) == 1)
                averageTheta = theta(i);
            else
                averageTheta = atan2(mean(sin(oldThetas{1, n}{1}(neighbours))), ...
                             mean(cos(oldThetas{1, n}{1}(neighbours))));
            end
            predictiveTheta = averageTheta + eta*r*deltaT;
            for j = 1:abs(h)
                predictiveTheta = predictiveTheta + averageTheta + eta*r*deltaT;
            end
            updatedTheta(i) = predictiveTheta + eta*r*deltaT;
        end
    else
        for i = 1:N
            r = rand - 0.5;
            neighbours = find(isNeighbour(i,:));
            if (length(neighbours) == 1)
                averageTheta = theta(i);
            else
                averageTheta = atan2(mean(sin(oldThetas{1, n-h}{1}(neighbours))), ...
                             mean(cos(oldThetas{1, n-h}{1}(neighbours))));
            end 
            updatedTheta(i) = averageTheta + eta*r*deltaT;
        end
    end
end