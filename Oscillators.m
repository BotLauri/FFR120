function cells = Oscillators(letterOfPattern, len, placement)
    if exist('placement', 'var')
        placement = placement - [1,1];
    end
    if ~exist('placement', 'var')
        placement = [0,0];
    end
    x = placement(1);
    y = placement(2);
    cells = zeros(len);
    if (letterOfPattern == 'a') % The blinker 1. 
        cells(1 + y, 1 + x) = 1;
        cells(2 + y, 1 + x) = 1;
        cells(3 + y, 1 + x) = 1;
    elseif (letterOfPattern == 'b') % The blinker 2.
        cells(1 + y, 1 + x) = 1;
        cells(1 + y, 2 + x) = 1;
        cells(1 + y, 3 + x) = 1;
    elseif (letterOfPattern == 'd') % The toad 1.
        cells(2 + y, 1 + x) = 1;
        cells(1 + y, 2 + x) = 1;
        cells(2 + y, 2 + x) = 1;
        cells(1 + y, 3 + x) = 1;
        cells(2 + y, 3 + x) = 1; 
        cells(1 + y, 4 + x) = 1;
    elseif (letterOfPattern == 'e' | letterOfPattern == 'g') % The toad 2 / The beacon 1.
        cells(1 + y, 1 + x) = 1;
        cells(2 + y, 1 + x) = 1;
        cells(1 + y, 2 + x) = 1;
        cells(3 + y, 4 + x) = 1;
        cells(4 + y, 3 + x) = 1;
        cells(4 + y, 4 + x) = 1; 
    elseif (letterOfPattern == 'h') % The beacon 2.
        cells(1 + y, 1 + x) = 1;
        cells(2 + y, 1 + x) = 1;
        cells(1 + y, 2 + x) = 1;
        cells(2 + y, 2 + x) = 1;
        cells(3 + y, 3 + x) = 1;
        cells(3 + y, 4 + x) = 1;
        cells(4 + y, 3 + x) = 1;
        cells(4 + y, 4 + x) = 1;
    end
end
