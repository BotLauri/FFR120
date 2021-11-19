function updatedTheta = UpdateOrientationWithDelay(theta, isNeighbour, eta, deltaT, h, oldThetas)
    S = 5; % Number of points for linear interpolation. 
    n = size(oldThetas, 2);
    if (h >= n)
        updatedTheta = theta;
        return
    end
    N = size(theta, 1);
    updatedTheta = zeros(N, 1);
    if (h < 0 && n > S) % Change this depending on what excerise (8.8 / 8.9).
        for i = 1:N
            neighbours = find(isNeighbour(i,:));
            averageTheta = zeros(S, 1);
            for s = 1:S
                if (length(neighbours) == 1)
                    averageTheta(s) = theta(i);
                else
                    averageTheta(s) = atan2(mean(sin(oldThetas{1, n-s}{1}(neighbours))), ...
                             mean(cos(oldThetas{1, n-s}{1}(neighbours))));
                end 
                averageTheta(s) = oldThetas{1, n-s}{1}(i);
            end
            poly = polyfit(1:S, averageTheta, 1);
            updatedTheta(i) = polyval(poly, S + abs(h));
        end
    elseif (h > 0)
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
    else
        updatedTheta = UpdateOrientation(theta, isNeighbour, eta, deltaT);
    end
end