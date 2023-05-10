function [blobs_] = seekCorners(blobs)
blobs_ = blobs(1,1:4);

UMIN = min(blobs(:).umin);
VMIN = min(blobs(:).vmin);
UMAX = max(blobs(:).umax);
VMAX = max(blobs(:).vmax);

k = 1;
for i = 1:length(blobs) % Checks which blob is:
    for j = 0:20
        if(blobs(1,i).umin == UMIN+j) % On the far left.
            blobs_(1,k) = blobs(1,i);
            k = k+1;
            break;
        end
        if(blobs(1,i).vmin == VMIN+j) % On the more above.
            blobs_(1,k) = blobs(1,i);
            k = k+1;
            break;
        end
        if(blobs(1,i).umax == UMAX-j) % On the far right.
            blobs_(1,k) = blobs(1,i);
            k = k+1;
            break;
        end
        if(blobs(1,i).vmax == VMAX-j) % On the more below. 
            blobs_(1,k) = blobs(1,i);
            k = k+1;
            break;
        end
    end
    if(k == 5) % After finding the 4 corners, return.
        break;
    end
end
end