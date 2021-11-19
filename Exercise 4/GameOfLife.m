%% Game of life. 
% Setup.
tic
clear all
len = 9;
maxNumberOfGenerations = 500;
isPeriodic = true; % Exercise 4.2c).
classification = zeros(1, 4); periodCounter = zeros(1, 5);
foundSpaceships = []; foundGliders = []; foundOscillators = [];

%deathOverPop = randi(8); deathUnderPop = randi(8); rebirth = randi(8);
%deathOverPop = 6; deathUnderPop = 2; rebirth = 5; % Fixed or oscillating stable pattern. 
%deathOverPop = 6; deathUnderPop = 4; rebirth = 5; % Extinction. 
% Complexity if deathOverPop >= deathUnderPop AND rebirth == 3.
% OR rebirth <= 2.

trials = 10^2; % Exercise 4.6b).
for t = 1:trials
cells = zeros(len);
isIdentical = false;

for i = floor(len/3):floor(2*len/3)
    for j = floor(len/3):floor(2*len/3)
        cells(i, j) = randi([0 1]); % Random initialization. 
    end
end

placement = [4, 4]; % Optional position parameter.
%cells = StillLife('a', len, placement);
%cells = Oscillators('b', len, placement); % Exercise 4.4.
%cells = Glider(3, len); % Exercise 4.5b).

% Exercise 4.3.
block = StillLife('a', len); 
beehive = StillLife('b', len); 
loaf = StillLife('c', len); 
boat = StillLife('d', len); 
tub = StillLife('e', len);
stillLife = [block, beehive, loaf, boat, tub];

blinker1 = Oscillators('a', len); blinker2 = Oscillators('b', len);
toad1 = Oscillators('d', len); toad2 = Oscillators('e', len);
beacon1 = Oscillators('g', len); beacon2 = Oscillators('h', len);
oscillators = [blinker1, blinker2, toad1, toad2, beacon1, beacon2];

northwest = Glider(1, len);
northeast = Glider(2, len);
southwest = Glider(3, len);
southeast = Glider(4, len);
gliders = [northwest, northeast, southwest, southeast];

llsp = Glider(5, len, placement);
%cells = llsp;

% Main loop. 
initialCells = cells;
previousCells = zeros(len);
newCells = zeros(len);
for n = 1:maxNumberOfGenerations
    %imagesc(cells) % Exercise 4.2b). 
    %pause(0.01)
    [~, cells] = GameOfLifeIteration(cells, len, isPeriodic);
    %[~, cells] = GameOfLifeIteration(cells, len, isPeriodic, deathOverPop, deathUnderPop, rebirth);
    if (isequal(cells, previousCells))
        break
    end
    if (n == maxNumberOfGenerations)
        break
    end
    previousCells = cells;
end

iterations = 0;
savedCells = cells;
while true
    iterations = iterations + 1;
    [savedCell, cells] = GameOfLifeIteration(cells, len, isPeriodic);
    %[savedCell, cells] = GameOfLifeIteration(cells, len, isPeriodic, deathOverPop, deathUnderPop, rebirth);
    savedCells = [savedCells, savedCell];
    for k = 1:iterations
        isFound = CheckConfiguration(cells, savedCells(:, ((k-1)*len+1):((k-1)*len+len)));
        if (isFound == true)
            break
        end
    end
    if (isFound == true)
        break
    end
end

% Classification: still life / dead, oscillators, gliders. 
% Iterations = 1 is always still life / dead.
% Iterations = len*period \leq 2*len for gliders.
% Therefore, iterations = 2 (or less than 2*len) is always at least one oscillator.
% The checks are needed as sometimes the maximum number of iterations is
% not enough for the game to converge.
% Because of these cases we need to check and throw away the ones which
% have not converged. 
% If many cases are thrown away -> Raise the maximum number of generations.
if (iterations == 1)
    classification(1) = classification(1) + 1;
    periodCounter(1) = periodCounter(1) + 1;
elseif (iterations >= 2*len)
    % Check if it really is a glider.
    period = iterations / len; foundCells = cells; isIdentical = false;   
    for m = 1:period
        [~, cells] = GameOfLifeIteration(cells, len, isPeriodic);
        %[~, cells] = GameOfLifeIteration(cells, len, isPeriodic, deathOverPop, deathUnderPop, rebirth);
    end
    for x = [-3, -2, -1, 0, 1, 2, 3]
        for y = [-3, -2, -1, 0, 1, 2, 3]
            if (~(x == 0 && y == 0))
                translatedCells = Translation(foundCells, x, y);
                isIdentical = CheckConfiguration(cells, translatedCells);
                if (isIdentical == true)
                    if (abs(x) == 1 && abs(y) == 1)
                        classification(3) = classification(3) + 1;
                        foundGliders = [foundGliders, cells];
                    elseif ((abs(x) == 2 && y == 0) || (x == 0 && abs(y) == 2))
                        classification(4) = classification(4) + 1;
                        foundSpaceships = [foundSpaceships, cells];
                    else
                        break
                    end
                    if (period > 4)
                        periodCounter(5) = periodCounter(5) + 1;
                    else
                        periodCounter(period) = periodCounter(period) + 1;
                    end
                    break
                end
            end
        end
    end
else
    % Check if it really is a oscillator.
    period = iterations; foundCells = cells; isIdentical = false;
    for m = 1:period
        [~, cells] = GameOfLifeIteration(cells, len, isPeriodic);
        %[~, cells] = GameOfLifeIteration(cells, len, isPeriodic, deathOverPop, deathUnderPop, rebirth);
    end
    isIdentical = CheckConfiguration(foundCells, cells);
    if (isIdentical == true)
        classification(2) = classification(2) + 1;
        foundOscillators = [foundOscillators, cells];
        if (period > 4)
            periodCounter(5) = periodCounter(5) + 1;
        else
            periodCounter(period) = periodCounter(period) + 1;
        end
    end
end

if (rem(t, 100) == 0)
    disp(t)
end

end
toc

%% Nice plot.
hold on
subplot(1,2,1)
imagesc(initialCells)
title('Initial Cells')
subplot(1,2,2)
imagesc(cells)
title('Cells')
hold off

%% Other plots.
imagesc(stillLife)
imagesc(oscillators)
imagesc(gliders)
imagesc(foundSpaceships(:,1:(0*len+len)))
imagesc(foundGliders(:,1:(1*len+len)))
imagesc(foundOscillators(:, 1:(2*len+len)))
