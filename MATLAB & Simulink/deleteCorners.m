function [im] = deleteCorners(im, corners)
UMIN = min(corners(:).umin);
VMIN = min(corners(:).vmin);
UMAX = max(corners(:).umax);
VMAX = max(corners(:).vmax);

for i = 1:length(corners) % Turn the corners values form 1 to 0.
    im(corners(i).vmin:corners(i).vmax,corners(i).umin:corners(i).umax) = 0;
end
im = im(VMIN:VMAX,UMIN:UMAX); % To reduce the image size to the workboard.
end

