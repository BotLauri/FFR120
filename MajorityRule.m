%% 4.8. Majority rule.
% Setup.
tic
clear all
len = 100;
p = 0.4;
isPeriodic = true;
cells = zeros(len);

trials = 10^1;
totalIterations = 0;
for t = 1:trials

for i = 1:len
    for j = 1:len
        r = rand();
        if (r < p)
            cells(i, j) = 1;
        else
            cells(i, j) = 0;
        end
    end
end

% Main loop. 
previousCells = zeros(len);
newCells = cells;
iterations = 0;
endProgram = false;
while (endProgram == false)
    for i = 1:len
        for j = 1:len
            neighboursVotes = AliveNeighbours(cells, i, j, len, isPeriodic);
            if (neighboursVotes == 0 | neighboursVotes == 1 | neighboursVotes == 2 | neighboursVotes == 3)
                newCells(i,j) = 0;
            elseif (neighboursVotes == 4)
            else
                newCells(i,j) = 1;
            end
        end
    end
    cells = newCells;
    iterations = iterations + 1;
    for k = 2:iterations
        if (isequal(cells, previousCells(:, ((k-1)*len+1):((k-1)*len+len))))
            endProgram = true;
        end
    end
    previousCells = [previousCells, cells];
end

totalIterations = totalIterations + iterations;
end 

%imagesc(cells)
avgIterations = totalIterations/trials
toc
