function yearsInPrison = PrisonersDilemma(m, n, N, R, S)
    if ~exist('S', 'var')
        S = 1.5;
    end
    % Initialization.
    T = 0; P = 1; 
    yearsInPrison = 0;
    playerA = [ones(1, m), zeros(1, N-m)];
    playerB = [ones(1, n), zeros(1, N-n)];
    
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
            yearsInPrison = yearsInPrison + P;
        elseif (strcmp(playerStrategies, '[0 1]')) % A defects, B cooperates.
            yearsInPrison = yearsInPrison + T;
        elseif (strcmp(playerStrategies, '[1 0]')) % A cooperates, B defects.
            yearsInPrison = yearsInPrison + S;
        else % Both cooperate.
            yearsInPrison = yearsInPrison + R;
        end
    end
end

