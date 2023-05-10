function [th] = seekBright(hist)
H = size(hist,1);
th = 0;

N = 100;
for i = H:-1:1
    if (hist(i,1) > N && th == 0)
        th = i/H; % Dummy. This just marks where the brightest part starts.
        continue;
    end
    if (hist(i,1) < N && th ~= 0)
        th = i/H; % Here the brightest part ends.
        break;
    end
end
end

% function [rang] = seekBright(hist)
% H = size(hist,1);
% rang = zeros(2, 1);
% 
% for i = H:-1:1
%     if (hist(i,1) > 100 && rang(2,1) == 0)
%         rang(2,1) = i/H;
%         continue;
%     end
%     if (hist(i,1) < 100 && rang(2,1) ~= 0)
%         rang(1,1) = i/H;
%         break;
%     end
% end
% end