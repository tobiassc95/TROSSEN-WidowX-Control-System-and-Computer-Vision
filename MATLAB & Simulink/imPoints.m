% This function sorts the points in the following order:
% 1--------2
% |        |
% |        |
% 4--------3

function [imP] = imPoints(imP)
if(size(imP,1) ~= 4)
    return;
end

for i = 1:3 % Sorts the points as explained.
    for j = (i+1):4
%         if(i==j) % No need to compare the same number.
%             continue;
%         end
        if(imP(i,1)>imP(j,1))
            imP([i j],:) = imP([j i],:);
        end
    end
end
%mP([2 3],:) = imP([3 2],:);

% Continue sorting the points as explained.
if(imP(2,2) < imP(1,2))
    imP([1 2],:) = imP([2 1],:);
end
imP([2 4],:) = imP([4 2],:);
if(imP(3,2) < imP(2,2))
    imP([2 3],:) = imP([3 2],:);
end
end

