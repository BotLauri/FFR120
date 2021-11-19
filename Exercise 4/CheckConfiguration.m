function isIdentical = CheckConfiguration(cells, configuration)
    isIdentical = true;
    if (size(cells, 1) ~= size(configuration, 1) || size(cells, 2) ~= size(configuration, 2))
        isIdentical = false;
        return
    end
    for i = 1:size(cells, 1)
        for j = 1:size(cells, 2)
            if (cells(i,j) ~= configuration(i,j))
                isIdentical = false;
                break
            end
        end
    end
end