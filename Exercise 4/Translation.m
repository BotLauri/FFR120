function newCells = Translation(cells, sx, sy)
    len = size(cells, 1);
    wid = size(cells, 2);
    newCells = zeros(len, wid);
    for i = 1:len
        for j = 1:wid
            if (cells(i,j) == 1)
                if (i + sy > len)
                    newY = i + sy - len;
                elseif (i + sy <= 0)
                    newY = i + sy + len;
                else
                    newY = i + sy;
                end
                if (j + sx > wid)
                    newX = j + sx - wid;
                elseif (j + sx <= 0)
                    newX = j + sx + wid;
                else
                    newX = j + sx;
                end
                newCells(newY, newX) = 1;
            end
        end
    end
end

