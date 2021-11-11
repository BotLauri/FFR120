function updatedTheta = UpdateOrientation(theta, isNeighbour, eta, deltaT)
    N = size(theta, 1);
    updatedTheta = zeros(N, 1);
    for i = 1:N
        updatedTheta(i) = sum(theta.*isNeighbour(:,i))/sum(isNeighbour(:,i)) + eta*randn*deltaT;
    end
end

