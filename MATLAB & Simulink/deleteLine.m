function [im] = deleteLine(im,blobs)
UMIN = min(blobs(:).umin);
VMIN = min(blobs(:).vmin);
UMAX = max(blobs(:).umax);
VMAX = max(blobs(:).vmax);

for i = VMIN:VMAX
    for j = UMIN:UMAX
        if (i<blobs(1,1).vmax && i>blobs(1,1).vmin)
            if (j<blobs(1,1).umax && j>blobs(1,1).umin)
                continue;
            end
        end
        if (i<blobs(1,2).vmax && i>blobs(1,2).vmin)
            if (j<blobs(1,2).umax && j>blobs(1,2).umin)
                continue;
            end
        end
        if (i<blobs(1,3).vmax && i>blobs(1,3).vmin)
            if (j<blobs(1,3).umax && j>blobs(1,3).umin)
                continue;
            end
        end
        if (i<blobs(1,4).vmax && i>blobs(1,4).vmin)
            if (j<blobs(1,4).umax && j>blobs(1,4).umin)
                continue;
            end
        end
        im(i,j) = 0;
    end
end
end

