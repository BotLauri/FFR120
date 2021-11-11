% Calculates the number of neighbours.
function amount = AliveNeighbours(cells, x, y, len, isPeriodic)
    if (isPeriodic == false)
        xNeighbours = max(1, x - 1):min(len, x + 1);
        yNeighbours = max(1, y - 1):min(len, y + 1);
    end
    if (isPeriodic == true)
        if (x == 1)
            if (y == 1)
                xNeighbours = [len, x, x + 1];
                yNeighbours = [len, y, y + 1];
            elseif (y == len)
                xNeighbours = [len, x, x + 1];
                yNeighbours = [y - 1, y, 1];
            else
                xNeighbours = [len, x, x + 1];
                yNeighbours = (y - 1):(y + 1);
            end
        elseif (x == len)
            if (y == 1)
                xNeighbours = [x - 1, x, 1];
                yNeighbours = [len, y, y + 1];
            elseif (y == len)
                xNeighbours = [x - 1, x, 1];
                yNeighbours = [y - 1, y, 1];
            else
                xNeighbours = [x - 1, x, 1];
                yNeighbours = (y - 1):(y + 1);
            end
        else
            if (y == 1)
                xNeighbours = (x - 1):(x + 1);
                yNeighbours = [len, y, y + 1];
            elseif (y == len)
                xNeighbours = (x - 1):(x + 1);
                yNeighbours = [y - 1, y, 1];
            else
                xNeighbours = (x - 1):(x + 1);
                yNeighbours = (y - 1):(y + 1);
            end
        end
    end
    amount = sum(sum(cells(xNeighbours, yNeighbours))) - cells(x, y);
end