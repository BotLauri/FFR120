function orientations = InitializeOrientations(N)
    orientations = zeros(N, 1);
    for i = 1:N
        orientations(i) = rand()*2*pi;
    end
end
