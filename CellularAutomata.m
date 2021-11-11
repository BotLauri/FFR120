%% Exercise 4.1. One-dimensional cellular automata. 
% b) Class 1, rule 160, all cells are 0/1. (uniform)
% Class 2, rule 178, stable/periodic structure. 
% Class 3, rule 182, chaotic pattern. (random)
% Class 4, rule 110, complicated structure. 
% Setup. 
clear all
wid = 150;
len = 100;
ruleNumberBase10 = input('Enter rule number: ');
if (ruleNumberBase10 < 0 || ruleNumberBase10 > 255)
    error('Rule must be between 0 and 255.')
end
ruleNumberBase2 = dec2bin(ruleNumberBase10, 8);
ruleNumberArray = zeros(1, 8);
ruleNumberArray(1) = str2double(ruleNumberBase2(1));
for i = 2:8
    ruleNumberArray(i) = str2double(ruleNumberBase2(i));
end

% Main loop.
cells = zeros(len, wid);
if (ruleNumberBase10 == 30 || ruleNumberBase10 == 90)
    cells(1, floor(wid/2) + 1) = 1; % One filled in the middle, Rule 30 and 90.
else
    for i = 1:wid
        cells(1,i) = randi([0 1]); % Random initialization, rule 184 and 110 (and others). 
    end
end
n = 2;
while (n <= len)
    for i = 1:wid
        if (i == 1)
            cells(n, 1) = UpdateRule(cells(n - 1, wid), cells(n - 1, 1), cells(n - 1, 2), ruleNumberArray);
        elseif (i == wid)
            cells(n, wid) = UpdateRule(cells(n - 1, wid - 1), cells(n - 1, wid), cells(n - 1, 1), ruleNumberArray);
        else
            cells(n, i) = UpdateRule(cells(n - 1, i - 1), cells(n - 1, i), cells(n - 1, i + 1), ruleNumberArray);
        end
    end
    n = n + 1;
end
imagesc(cells)

function output = UpdateRule(a, b, c, ruleNumberArray)
    if (a == 1)
        if (b == 1)
            if (c == 1)
                output = ruleNumberArray(1);
            else
                output = ruleNumberArray(2);
            end
        else
            if (c == 1)
                output = ruleNumberArray(3);
            else
                output = ruleNumberArray(4);
            end
        end
    else
        if (b == 1)
            if (c == 1)
                output = ruleNumberArray(5);
            else
                output = ruleNumberArray(6);
            end
        else
            if (c == 1)
                output = ruleNumberArray(7);
            else
                output = ruleNumberArray(8);
            end
        end
    end
end
