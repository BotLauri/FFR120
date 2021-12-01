%% Prisoner's dilemma with multiple rounds. 
clear all
tic
N = 10; T = 0; R = 0.5; P = 1; S = 1.5;

yearsInPrison = zeros(1, N+1);
for t = 0:N
    % Initialization. 
    n = t; m = 6;
    playerA = [ones(1, t), zeros(1, N-t)];
    playerB = [ones(1, m), zeros(1, N-m)];
    
    % Game loop.
    for i = 1:N
        indexA = find(playerA == 0, 1, 'first');
        if (isempty(indexA))
            indexA = N + 1;
        end
        indexB = find(playerB == 0, 1, 'first');
        if (isempty(indexB))
            indexB = N + 1;
        end
        if (indexA - 1 > indexB)
            playerA(indexA - 1) = 0;
        end
        if (indexB - 1 > indexA)
            playerB(indexB - 1) = 0;
        end
    end
    
    % Calculate years in prison. 
    for i = 1:N
        playerStrategies = mat2str([playerA(i), playerB(i)]);
        if (strcmp(playerStrategies, '[0 0]')) % Both defect.
            yearsInPrison(t+1) = yearsInPrison(t+1) + P;
        elseif (strcmp(playerStrategies, '[0 1]')) % A defects, B cooperates.
            yearsInPrison(t+1) = yearsInPrison(t+1) + T;
        elseif (strcmp(playerStrategies, '[1 0]')) % A cooperates, B defects.
            yearsInPrison(t+1) = yearsInPrison(t+1) + S;
        else % Both cooperate.
            yearsInPrison(t+1) = yearsInPrison(t+1) + R;
        end
    end
end

hold on
x = [6 6]; y = [6 10.5];
line(x, y,'Color', 'black', 'LineStyle', '--')
scatter(0:N, yearsInPrison, 'filled')
ylim([6, 9.5])
xlabel('n')
ylabel('Years in prison')
hold off

toc
