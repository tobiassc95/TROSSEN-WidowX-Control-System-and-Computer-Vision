function [linP] = seekLine(im,rho,theta)
linP = zeros(2,2);
V = size(im,1);
U = size(im,2);

if abs(theta) <= pi/4
    v = zeros(1,U);
    for i = 1:U
        v(i) = round(-i*tan(theta)+(rho/cos(theta))); % Line equation from hough parameters.
    end
    v = v+1; % Since array index starts with 1 in Matlab.

    for i = 1:U % Seachs the points of the line.
        if (v(i) > 0 && v(i) <= V)
            if(im(v(i),i)==1 && linP(1,1)==0)
                linP(1,:) = [i*200/U v(i)*150/V];
            end
            if(im(v(i),i)==0 && linP(1,1)~=0)
                linP(2,:) = [(i-1)*200/U v(i-1)*150/V];
                break;
            end
        end
    end
else %if abs(theta) > pi/4
    u = zeros(1,V);
    for i = 1:V
        u(i) = round(-i/tan(theta)+(rho/sin(theta))); % Line equation from hough parameters.
    end
    u = u+1; % Since array index starts with 1 in Matlab.

    for i = 1:V % Seachs the points of the line.
        if (u(i) > 0 && u(i) <= U)
            if(im(i,u(i))==1 && linP(1,1)==0)
                linP(1,:) = [u(i)*200/U i*150/V];
            end
            if(im(i,u(i))==0 && linP(1,1)~=0)
                linP(2,:) = [u(i-1)*200/U (i-1)*150/V];
                break;
            end
        end
    end
end

% Sorts the points from the closest (to the manipulator) to the farest.
if(linP(1,2) < linP(2,2))
    linP([1 2],:) = linP([2 1],:);
% elseif(linP(1,2) == linP(2,2))
%     if(linP(1,1) > linP(2,1))
%         linP([1 2],:) = linP([2 1],:);
%     end
end
end

