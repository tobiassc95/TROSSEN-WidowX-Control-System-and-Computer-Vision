function [im] = deleteBackground(im)
V = size(im,1);
U = size(im,2);

for i = 1:V % Scans the image and turn the background values form 1 to 0.
    if(im(i,1) == 1)
        for j = 1:U
            if (im(i,j) == 0)
                break;
            end
            im(i,j) = 0;
        end
    end
    if(im(i,U) == 1)
        for j = U:-1:1
            if (im(i,j) == 0)
                break;
            end
            im(i,j) = 0;
        end
    end
end
for i = 1:U % Scans the image and turn the background values form 1 to 0.
    if(im(1,i) == 1)
        for j = 1:V
            if (im(j,i) == 0)
                break;
            end
            im(j,i) = 0;
        end
    end
    if(im(V,i) == 1)
        for j = V:-1:1
            if (im(j,i) == 0)
                break;
            end
            im(j,i) = 0;
        end
    end
end
end

% function [im] = deleteBackground(im)
% V = size(im,1);
% U = size(im,2);
% 
% for i = 1:V % Scans the image and turn the background values form 1 to 0.
%     if(im(i,1) == 1)
%         for j = 1:U
%             if (im(i,j)==0 && j~=U)
%                 if(im(i,j+1) == 0)
%                     break;
%                 end
%             end
%             im(i,j) = 0;
%         end
%     end
%     if(im(i,U) == 1)
%         for j = U:-1:1
%             if (im(i,j)==0 && j~=1)
%                 if(im(i,j-1) == 0)
%                     break;
%                 end
%             end
%             im(i,j) = 0;
%         end
%     end
% end
% for i = 1:U % Scans the image and turn the background values form 1 to 0.
%     if(im(1,i) == 1)
%         for j = 1:V
%             if (im(j,i)==0 && j~=V)
%                 if(im(j+1,i) == 0)
%                     break;
%                 end
%             end
%             im(j,i) = 0;
%         end
%     end
%     if(im(V,i) == 1)
%         for j = V:-1:1
%             if (im(j,i)==0 && j~=1)
%                 if(im(j-1,i) == 0)
%                     break;
%                 end
%             end
%             im(j,i) = 0;
%         end
%     end
% end
% end