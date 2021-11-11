function [previousCells, cells] = GameOfLifeIteration(cells, len, isPeriodic, deathOverPop, deathUnderPop, rebirth)
    if ~exist('deathOverPop', 'var')
        deathOverPop = 3;
    end
    if ~exist('deathUnderPop', 'var')
        deathUnderPop = 2;
    end
    if ~exist('rebirth', 'var')
        rebirth = 3;
    end
    previousCells = cells;
    newCells = zeros(len);
    for i = 1:len
        for j = 1:len
            aliveNeighbours = AliveNeighbours(cells, i, j, len, isPeriodic);
            if (aliveNeighbours > deathOverPop && cells(i,j) == 1)
                newCells(i,j) = 0;
            elseif (aliveNeighbours < deathUnderPop && cells(i,j) == 1)
                newCells(i,j) = 0;
            elseif (aliveNeighbours == rebirth && cells(i,j) == 0)
                newCells(i,j) = 1;
            else
                newCells(i,j) = cells(i,j);
            end
        end
    end
    cells = newCells;
end
