function [im] = deleteIslands(im)
V = size(im,1);
U = size(im,2);

W = 15;
for i = 1:V-W
    for j = 1:U-W
        edge = false;
        for k = i:i+W
            if(im(k,j) == 1)
                edge = true;
                break;
            end
        end
        for k = i:i+W
            if(im(k,j+W) == 1)
                edge = true;
                break;
            end
        end
        for k = j:j+W
            if(im(i,k) == 1)
                edge = true;
                break;
            end
        end
        for k = j:j+W
            if(im(i+W,k) == 1)
                edge = true;
                break;
            end
        end
        if(~edge)
            im(i:i+W,j:j+W) = 0;
        end
    end
end
end

