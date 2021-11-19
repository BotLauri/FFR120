function updatedTheta = UpdateOrientationNearestNeighbours(theta, r, k, eta, deltaT, vision)
    N = size(theta, 1);
    updatedTheta = zeros(N, 1);
    for i = 1:N
        distances = zeros(N, 1);
        for n = 1:N
            distances(n) = Inf;
        end
        for j = 1:N
            angle = acosd(dot(r(i, :), r(j, :))/(norm(r(i, :))*norm(r(j, :))));
            if (angle < vision)
                distances(j) = norm(r(i, :) - r(j, :));
            end
        end
        [~, indexes] = sort(distances);
        nearestNeighbours = indexes(1:k);
        averageTheta = atan2(mean(sin(theta(nearestNeighbours))), mean(cos(theta(nearestNeighbours))));
        updatedTheta(i) = averageTheta + eta*(rand - 0.5)*deltaT;
    end
end