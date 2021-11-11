function particlePositions = InitializePositions(N, L)
    particlePositions = zeros(N, 2);
    for i = 1:N
        r = rand(); 
        q = rand();
        x = L*r - L/2;
        y = L*q - L/2;
        particlePositions(i, :) = [x, y];
    end
end

