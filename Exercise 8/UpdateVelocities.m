function velocities = UpdateVelocities(v, theta)
    N = size(theta, 1);
    velocities = zeros(N, 2);
    for i = 1:N
        velocities(i, 1) = v*(cos(theta(i)));
        velocities(i, 2) = v*(sin(theta(i)));
    end
end
