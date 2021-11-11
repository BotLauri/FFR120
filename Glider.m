function cells = Glider(direction, len, placement)
    if exist('placement', 'var')
        placement = placement - [1,1];
    end
    if ~exist('placement', 'var')
        placement = [0,0];
    end
    x = placement(1);
    y = placement(2);
    cells = zeros(len);
    if (direction == 1) % Northwest. 
        cells(1 + y, 1 + x) = 1;
        cells(2 + y, 1 + x) = 1;
        cells(3 + y, 2 + x) = 1;
        cells(1 + y, 2 + x) = 1;
        cells(1 + y, 3 + x) = 1;
    elseif (direction == 2) % Northeast. 
        cells(1 + y, 1 + x) = 1;
        cells(2 + y, 3 + x) = 1;
        cells(3 + y, 2 + x) = 1;
        cells(1 + y, 2 + x) = 1;
        cells(1 + y, 3 + x) = 1;
    elseif (direction == 3) % Southwest. 
        cells(2 + y, 1 + x) = 1;
        cells(1 + y, 2 + x) = 1;
        cells(3 + y, 1 + x) = 1;
        cells(3 + y, 2 + x) = 1;
        cells(3 + y, 3 + x) = 1;
    elseif (direction == 4) % Southeast. 
        cells(2 + y, 3 + x) = 1;
        cells(1 + y, 2 + x) = 1;
        cells(3 + y, 1 + x) = 1;
        cells(3 + y, 2 + x) = 1;
        cells(3 + y, 3 + x) = 1;
    elseif (direction == 5) % Lightweight spaceship.
        cells(2 + y, 3 + x) = 1;
        cells(2 + y, 6 + x) = 1;
        cells(3 + y, 2 + x) = 1;
        cells(4 + y, 2 + x) = 1;
        cells(4 + y, 6 + x) = 1;
        cells(5 + y, 2 + x) = 1;
        cells(5 + y, 3 + x) = 1;
        cells(5 + y, 4 + x) = 1;
        cells(5 + y, 5 + x) = 1;
    end        
end
