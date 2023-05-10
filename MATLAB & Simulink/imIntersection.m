function [intP] = imIntersection(imL1,imL2)
    intP = zeros(1,2);
    
    imL2 = 2*imL2;
    imL = imL1+imL2;
    [v,u] = find(imL);
    
    points = size(v);
    %count = 0;
    %intersection = 0;
    for i = 1:points
        count = 0;
        if(v(i)>1 && v(i)<max(v))
            if(u(i)>1 && u(i)<max(u))
                for j = v(i)-1:v(i)+1
                    for k = u(i)-1:u(i)+1
                        if(imL(j,k)==0)
                            continue;
                        elseif(imL(j,k)==3)
                            intP(1) = j;
                            intP(2) = k;
                            return;
                        else
                            if(imL(j,k)~=imL(v(i),u(i)))
                                count = count+1;
                                if(count > 1)
                                    %imL(v(i),u(i)) = 3;
                                    intP = [v(i) u(i)];
                                    return;
                                end
                            end
                        end
                    end
                end
            end
        end
	end
end
