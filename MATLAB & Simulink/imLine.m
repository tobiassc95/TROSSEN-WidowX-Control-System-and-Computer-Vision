function [I_] = imLine(I,rho,theta)
V = size(I,1);
U = size(I,2);
I_ = zeros(V,U);

v = zeros(1,U);
for i = 1:U
    v(i) = round(-i*tan(theta)+(rho/cos(theta)));
end

v = v+1;
for i = 1:U
    if (v(i) > 0 && v(i) <= V)
        I_(v(i),i) = 1;
    end
end

% % Mirror. To fill more pixels in the line.
% 
% u = zeros(1,V);
% for i = 1:V
% 	u(i) = round(-(i/tan(theta))+(rho/sin(theta)));
% end
% 
% u = u+1;
% for i = 1:V
%     if (u(i) > 0 && u(i) <= U)
%         I_(i,u(i)) = 1;
%     end
% end
end