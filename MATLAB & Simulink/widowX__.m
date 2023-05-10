clear
close all
clc

%% Information.
workboardX = 0.15; % Width.
workboardY = 0.2; % Length.
workboardZ = 0.03; % Height.
error = 0.03; % Error in the 'x' axis due to the model.
Xoffset = 0.12; % Distance between the origin and the edge of the workboard.
toolLength = 0.03; %Length of the tool.

%% Image noise filtering.
%im = iread('IMG_7702.jpg');
%im = iread('IMG_7703.jpg');
%im = iread('IMG_7786.jpg'); % GOOD!
%im = iread('IMG_7710.jpg'); % PROBLEMS WITH 7781,7709,(7710),(7711),7712.

% VERTICAL
% im = iread('IMG_0594.jpg'); % GOOD!
% im = iread('IMG_0600.jpg'); % GOOD!

% HORIZONTAL
% im = iread('IMG_0603.jpg'); % GOOD!
% im = iread('IMG_0602.jpg'); % BAD!
% im = iread('IMG_0601.jpg'); % BAD/BAD!

% LONG DIAGONAL
% im = iread('IMG_0606.jpg'); % GOOD!
% im = iread('IMG_0607.jpg'); % BAD/GOOD!

% SHORT DIAGONAL
im = iread('IMG_0609.jpg'); % GOOD!
% im = iread('IMG_0610.jpg'); % GOOD!

im = idouble(im);
about im;
im = iconvolve(im,kgauss(1)); % To erase noise in the image.
figure(1);
idisp(im);

%% Image in grayscale.
im = imono(im);
figure(2);
idisp(im);

%% Image brightness filtering.
th = seekBright(ihist(im)); % To separate the brightess part (paper) with the rest.
% im_ = im;
im = im < th(1);
figure(3);
idisp(im);

%% Background removing
im = deleteBackground(im); % Now the background will be black.
% figure(4);
% idisp(im);

%% Noise left filtering
im = deleteIslands(im);

% S = ones(2);
% S = kcircle(1);
% %imG = irank(imG, 3, S);
% %im = iclose(im,S);
% im = iopen(im,S);
figure(4);
idisp(im);

%% BLOBS
blobs = iblobs(im); % To find figures withing the brightest part (paper) of the image.
blobs = blobs(blobs.touch == 0); % Remove the blob containing the backgroung.
%blobs = blobs(blobs.parent ~= 0);
corners = seekCorners(blobs); % To store the 4 corners and remove the other blobs.

hold on;
for i=1:length(corners) % To see the corners and their center of mass (COM).
    plot_box([corners(i).umin corners(i).umax; corners(i).vmin corners(i).vmax], 'r');
    plot_box([corners(i).uc-1 corners(i).uc+1; corners(i).vc-1 corners(i).vc+1], 'g');
end

%% WRAPING
hold off;
V = size(im,1);
U = size(im,2);

imP = zeros(4,2);
for i = 1:length(corners) % Stores the COM of each blob.
    imP(i,:) = [corners(i).uc corners(i).vc];
end
imP = imPoints(imP);
imP = imP'; 
imP_ = [1 U U 1;1 1 V V]; % Destination points to do the warpng. 

matH = homography(imP,imP_); % To generate the homography matrix.
imW = homwarp(matH,im,'full'); % To change the perspective of the image.
imW = imW>0.5;
figure(5);
idisp(imW)

%% Remove corners.
blobs = iblobs(imW);
blobs = blobs(blobs.touch == 0); % Remove the blob containing the backgroung.
%blobs = blobs(blobs.parent ~= 0);
corners = seekCorners(blobs); % To store the 4 corners and remove the other blobs.

hold on;
for i=1:length(corners) % To see the corners and their center of mass (COM).
    plot_box([corners(i).umin corners(i).umax; corners(i).vmin corners(i).vmax], 'r');
    plot_box([corners(i).uc-1 corners(i).uc+1; corners(i).vc-1 corners(i).vc+1], 'g');
end

hold off;
im = deleteCorners(imW,corners); % To remove the corners from the image.
about im;
figure(9);
idisp(im);

%% Line information
imH = Hough(im,'nbins',[800 802]); % To find the line and its equation.
%imH = Hough(im,'nbins',[V U]);
%imH = Hough(im);
%imH.edgeThresh = 0;
imH.houghThresh = 0.5;
imH.suppress = 80;
imH.plot;
% figure(10)
% imH.lines;
lnH = imH.lines;

linP = seekLine(im,lnH(1).rho,lnH(1).theta)/1000;  %linP = [yi xi; yf xf];
linP = [workboardY/2-linP(1,1) Xoffset+(workboardX-linP(1,2));
        workboardY/2-linP(2,1) Xoffset+(workboardX-linP(2,2))]
linP_ = linP + [0 -error;
                0 -error];

% linP = [0.05 workboardX; 0.05 0.25];

%% Robot (The whole robot).
N = 6; %Links.
L = [0.13 0.144 0.056 0.144 0.144/2 0.144/2+toolLength]; %Links length.
%m = [1 1]; %Links weight.

DH = struct('alpha', cell(1,N), 'a', cell(1,N), 'd', cell(1,N)); % Structure.
DH(1).alpha = 0;        DH(1).a = 0;    DH(1).d = L(1); %DH(1).type = 'R';
DH(2).alpha = pi/2;     DH(2).a = 0;    DH(2).d = 0; %DH(2).type = 'R';
DH(3).alpha = -pi/2;    DH(3).a = 0;    DH(3).d = L(2); %DH(1).type = 'R';
DH(4).alpha = pi/2;   DH(4).a = L(3); DH(4).d = 0; %DH(2).type = 'R';
DH(5).alpha = 0;      DH(5).a = L(4); DH(5).d = 0; %DH(1).type = 'R';
DH(6).alpha = pi/2;   DH(6).a = 0;    DH(6).d = L(5); %DH(2).type = 'R';
%DH(7).alpha = 0;       DH(7).a = 0;    DH(7).d = L(6); %DH(2).type = 'R';

links = {0};
% for i = 1:N
%   links{i} = Link('alpha', DH(i).alpha, 'a', DH(i).a, 'd', DH(i).d, ...
%       'modified', 'm', m(i), 'r', [L(i),0,0], 'G', 1, 'Tc', 1);%'Tc', [1 1]); % Links vector.
% 	%links{i}.dyn();
% end
for i = 1:N
    links{i} = Link('alpha', DH(i).alpha, 'a', DH(i).a, 'd', DH(i).d, ...
        'modified');%, 'G', 1); % Links vector.
    %links{i}.dyn();
end
size(links)

tool = transl([0, 0, L(6)]); % Tool offset.

widowXbot = SerialLink([links{:}], 'tool', tool, 'name', 'widowXbot');
%widowXbot = SerialLink([links{:}], 'name', 'widowXbot');
%widowXbot.dyn();
%RRbotD = RRbot.perturb(0.8); %Disturbance.

q0=[0 0 0 0 0 0]; %Inintial angle.
%figure();
%widowXbot.teach(q0); %q3 is not actually present.
% hold on;
% patch( [2 0 0 2] , [0 2 2 0], [5 5 -5 -5], 'white')
% %patch( [1 -1 -1 1] , [0 0 0 0], [1 1 -1 -1], [1 1 -1 -1])

%% Robot_ (A simpler model of the robot).
N_ = 3;
L1_ = sqrt(L(2)^2+L(3)^2);
L_ = [L(1) L1_ L(4)]; %Links length.

DH_ = struct('alpha', cell(1,N_), 'a', cell(1,N_), 'd', cell(1,N_)); % Structure.
DH_(1).alpha = 0;        DH_(1).a = 0;    DH_(1).d = L_(1); %DH_(1).type = 'R';
DH_(2).alpha = pi/2;     DH_(2).a = 0;    DH_(2).d = 0; %DH_(2).type = 'R';
%DH_(3).alpha = -pi/2;    DH_(3).a = 0;    DH_(3).d = L(2); %DH_(1).type = 'R';
DH_(3).alpha = 0;   DH_(3).a = L_(2); DH_(3).d = 0; %DH_(2).type = 'R';
% DH_(4).alpha = 0;      DH_(4).a = L(4); DH_(4).d = 0; %DH(1).type = 'R';
% DH_(5).alpha = pi/2;   DH_(5).a = 0;    DH_(5).d = L(5); %DH(2).type = 'R';
% %DH_(6).alpha = 0;       DH_(6).a = 0;    DH_(6).d = L(6); %DH(2).type = 'R';

links_ = {0};
for i = 1:N_
    links_{i} = Link('alpha', DH_(i).alpha, 'a', DH_(i).a, 'd', DH_(i).d, ...
        'modified');%, 'G', 1); % Links vector.
    %links{i}.dyn();
end
size(links_)

%tool_ = transl([L(4), 0, 0]); % Tool offset.
tool_ = transl([L_(3), 0, 0]); % Tool offset.

widowXbot_ = SerialLink([links_{:}], 'tool', tool_, 'name', 'widowXbot_');
%widowXbot_.dyn();

q0_ = [0 pi/2-atan2(L(3),L(2)) -(pi/2-atan2(L(3),L(2)))];
%figure();
%widowXbot_.teach(q0_); %q3 is not actually present.

%% Robot__
N__ = 5;
L__ = [L(1) L1_ L(4) L(5) L(6)]; %Links length.

% A simpler model of the robot.
DH__ = struct('alpha', cell(1,N_), 'a', cell(1,N_), 'd', cell(1,N_)); % Structure.
DH__(1).alpha = 0;        DH__(1).a = 0;    DH__(1).d = L__(1); %DH_(1).type = 'R';
DH__(2).alpha = pi/2;     DH__(2).a = 0;    DH__(2).d = 0; %DH_(2).type = 'R';
%DH__(3).alpha = -pi/2;    DH__(3).a = 0;    DH__(3).d = L(2); %DH_(1).type = 'R';
DH__(3).alpha = 0;         DH__(3).a = L__(2); DH__(3).d = 0; %DH__(2).type = 'R';
DH__(4).alpha = 0;      DH__(4).a = L__(3); DH__(4).d = 0; %DH(1).type = 'R';
DH__(5).alpha = pi/2;   DH__(5).a = 0;    DH__(5).d = L__(4); %DH(2).type = 'R';
% %DH__(6).alpha = 0;       DH__(6).a = 0;    DH__(6).d = L(6); %DH(2).type = 'R';

links__ = {0};
for i = 1:N__
    links__{i} = Link('alpha', DH__(i).alpha, 'a', DH__(i).a, 'd', DH__(i).d, ...
        'modified');%, 'G', 1); % Links vector.
    %links{i}.dyn();
end
size(links__)

tool__ = transl([0, 0, L__(5)]); % Tool offset.

widowXbot__ = SerialLink([links__{:}], 'tool', tool__, 'name', 'widowXbot__');
%widowXbot__.dyn();

q0__ = [0 pi/2-atan2(L(3),L(2)) -(pi/2-atan2(L(3),L(2))) 0 0];
% figure();
%widowXbot__.teach(q0__); %q3 is not actually present.

%% Time Animation.
tmax = 1; %time max.
ts = 0.05; %Time step.
steps = 10;

%% Cartesian spaces (scene 1).
Ti = transl(L(3)+L(4), 0, L(1)+L(2)-L(5)-L(6)+L(1));%(0.25, 0.05, L(1));%L(5)+L(6)+workboardZ);
Tf = transl(linP_(1,2), linP_(1,1), L(1)+workboardZ);%L(5)+L(6)+workboardZ);
T = ctraj(Ti, Tf, steps);

Ti = transl(L(3)+L(4), 0, L(1)+L(2)-L(5)-L(6));%(0.25, 0.05, 0);%-L(1)+L(5)+L(6)+workboardZ);
Tf = transl(linP_(1,2), linP_(1,1), workboardZ);%-L(1)+L(5)+L(6)+workboardZ);
T_ = ctraj(Ti, Tf, steps);
x = [T_(1,4,:);T_(2,4,:);T_(3,4,:)]; % Gets the the position.

%% Angular Coordinates (scene 1).
q1 = atan2(x(2,:),x(1,:));
q3 = asin((x(1,:).^2+x(2,:).^2+x(3,:).^2-L_(2)^2-L_(3)^2)/(2*L_(2)*L_(3)));
q2_ = atan2(L(3),L(2));
q2 = asin((x(3,:)*L(4).*cos(q3)-(L1_+L(4)*sin(q3)).*sqrt(x(1,:).^2+x(2,:).^2))./(L1_^2+2*L1_*L(4)*sin(q3)+L(4)^2))+q2_;
qs1 = [q1; q2; zeros(1, steps); q3; -q2-q3; zeros(1, steps)];
size(qs1)

%q_ = widowXbot_.ikine(T_, 'mask',[1 1 1 0 0 0])';
q_ = widowXbot_.ikine(T,'q0', q0_, 'mask',[1 1 1 0 0 0])';
qs1_ = [q_(1,:); q_(2,:)-pi/2+q2_; zeros(1,steps); q_(3,:)+pi/2; -q_(2,:)-q2_-q_(3,:); zeros(1,steps)];
size(qs1_)

qs1-qs1_

%% Cartesian spaces (scene 2).
Ti = transl(linP_(1,2), linP_(1,1), L(1)+workboardZ);%(0.25, 0.05, L(1));%L(5)+L(6)+workboardZ);
Tf = transl(linP_(2,2), linP_(2,1), L(1)+workboardZ);%L(5)+L(6)+workboardZ);
T = ctraj(Ti, Tf, steps);

Ti = transl(linP_(1,2), linP_(1,1), workboardZ);%(0.25, 0.05, 0);%-L(1)+L(5)+L(6)+workboardZ);
Tf = transl(linP_(2,2), linP_(2,1), workboardZ);%-L(1)+L(5)+L(6)+workboardZ);
T_ = ctraj(Ti, Tf, steps);
x = [T_(1,4,:);T_(2,4,:);T_(3,4,:)]; % Gets the the position.

%% Angular Coordinates (scene 2).
q1 = atan2(x(2,:),x(1,:));
q3 = asin((x(1,:).^2+x(2,:).^2+x(3,:).^2-L_(2)^2-L_(3)^2)/(2*L_(2)*L_(3)));
q2_ = atan2(L(3),L(2));
q2 = asin((x(3,:)*L(4).*cos(q3)-(L1_+L(4)*sin(q3)).*sqrt(x(1,:).^2+x(2,:).^2))./(L1_^2+2*L1_*L(4)*sin(q3)+L(4)^2))+q2_;
qs2 = [q1; q2; zeros(1, steps); q3; -q2-q3; zeros(1, steps)];
size(qs2)

%q_ = widowXbot_.ikine(T_, 'mask',[1 1 1 0 0 0])';
q_ = widowXbot_.ikine(T,'q0', q0_, 'mask',[1 1 1 0 0 0])';
qs2_ = [q_(1,:); q_(2,:)-pi/2+q2_; zeros(1,steps); q_(3,:)+pi/2; -q_(2,:)-q2_-q_(3,:); zeros(1,steps)];
size(qs2_)

qs2-qs2_

%% Cartesian spaces (scene 3).
Ti = transl(linP_(2,2), linP_(2,1), L(1)+workboardZ);%(0.25, 0.05, L(1));%L(5)+L(6)+workboardZ);
Tf = transl(L(3)+L(4), 0, L(1)+L(2)-L(5)-L(6)+L(1));%L(5)+L(6)+workboardZ);
T = ctraj(Ti, Tf, steps);

Ti = transl(linP_(2,2), linP_(2,1), workboardZ);%(0.25, 0.05, 0);%-L(1)+L(5)+L(6)+workboardZ);
Tf = transl(L(3)+L(4), 0, L(1)+L(2)-L(5)-L(6));%-L(1)+L(5)+L(6)+workboardZ);
T_ = ctraj(Ti, Tf, steps);
x = [T_(1,4,:);T_(2,4,:);T_(3,4,:)]; % Gets the the position.

%% Angular Coordinates (scene 3).
q1 = atan2(x(2,:),x(1,:));
q3 = asin((x(1,:).^2+x(2,:).^2+x(3,:).^2-L_(2)^2-L_(3)^2)/(2*L_(2)*L_(3)));
q2_ = atan2(L(3),L(2));
q2 = asin((x(3,:)*L(4).*cos(q3)-(L1_+L(4)*sin(q3)).*sqrt(x(1,:).^2+x(2,:).^2))./(L1_^2+2*L1_*L(4)*sin(q3)+L(4)^2))+q2_;
qs3 = [q1; q2; zeros(1, steps); q3; -q2-q3; zeros(1, steps)];
size(qs3)

%q_ = widowXbot_.ikine(T_, 'mask',[1 1 1 0 0 0])';
q_ = widowXbot_.ikine(T,'q0', q0_, 'mask',[1 1 1 0 0 0])';
qs3_ = [q_(1,:); q_(2,:)-pi/2+q2_; zeros(1,steps); q_(3,:)+pi/2; -q_(2,:)-q2_-q_(3,:); zeros(1,steps)];
size(qs3_)

qs3-qs3_

%% Workspace
DH_pw(1) = Link([0 0 L1_ 0]);
DH_pw(2) = Link([0 0 L(4) 0]);
q2_pw = (-pi/2+pi/6:0.05:pi/6);
q3_pw = (-pi/4:-0.05:-pi/2-pi/4);
q_pw = {q2_pw, q3_pw};

figure(10);
plotworkspace(DH_pw, q_pw);
hold on;
rectangle('Position',[0.1,0,0.15,0.03],'FaceColor',[1 1 1], 'EdgeColor','black');
yline(0,'black');
% line([-0.1;0.01;0.3],[-0.01;0.001;0.01],'black');
% rectangle('Position',[-0.1,-0.001,0.3,0.001],'FaceColor',[0 0 0]);
% x_pw = (-0.1:0.01:0.3);
% y_pw = zeros(size(x_pw));
% plot(x_pw,y_pw,'black');

%% Animation
hold off;
figure(11);
widowXbot.plot(q0);
hold on;
patch([2 -2 -2 2], [2 2 -2 -2], [0 0 0 0], "black");
patch([Xoffset Xoffset Xoffset+workboardX Xoffset+workboardX], [-workboardY/2 workboardY/2 workboardY/2 -workboardY/2], [workboardZ workboardZ workboardZ workboardZ], "white");
patch([Xoffset Xoffset Xoffset+workboardX Xoffset+workboardX], [-workboardY/2 workboardY/2 workboardY/2 -workboardY/2], [0 0 0 0], "white");
patch([Xoffset Xoffset Xoffset Xoffset], [-workboardY/2 workboardY/2 workboardY/2 -workboardY/2], [0 0 workboardZ workboardZ], "white");
patch([Xoffset+workboardX Xoffset+workboardX Xoffset+workboardX Xoffset+workboardX], [-workboardY/2 workboardY/2 workboardY/2 -workboardY/2], [0 0 workboardZ workboardZ], "white");
patch([Xoffset Xoffset+workboardX Xoffset+workboardX Xoffset], [-workboardY/2 -workboardY/2 -workboardY/2 -workboardY/2], [0 0 workboardZ workboardZ], "white");
patch([Xoffset Xoffset+workboardX Xoffset+workboardX Xoffset], [workboardY/2 workboardY/2 workboardY/2 workboardY/2], [0 0 workboardZ workboardZ], "white");
patch([linP(1,2)-0.003 linP(1,2)+0.003 linP(2,2)+0.003 linP(2,2)-0.003], [linP(1,1)-0.003 linP(1,1)+0.003 linP(2,1)+0.003 linP(2,1)-0.003], [workboardZ workboardZ workboardZ workboardZ], "red");
patch([linP(1,2)-0.003 linP(1,2)+0.003 linP(2,2)+0.003 linP(2,2)-0.003], [linP(1,1)-0.003 linP(1,1)+0.003 linP(2,1)+0.003 linP(2,1)-0.003], [0 0 0 0], "red");
%linP = [yi xi; yf xf];
while 1
    pause(1);
    
    widowXbot.plot(qs1');
    %widowXbot.plot(qs1_');
    %widowXbot__.plot([q1; q2+pi/2; q3-pi/2; -q2-q3; zeros(1,steps)]');
    
    widowXbot.plot(qs2');
    %widowXbot.plot(qs2_');
    %widowXbot__.plot([q1; q2+pi/2; q3-pi/2; -q2-q3; zeros(1,steps)]');
    
    widowXbot.plot(qs3');
    %widowXbot.plot(qs3_');
    %widowXbot__.plot([q1; q2+pi/2; q3-pi/2; -q2-q3; zeros(1,steps)]');
end

%% Other stuff.
% invJ = zeros(3,3,steps);
% for i = 1:steps
%     invJ(:,:,i) = [L1_*sin(q(1,i)).*sin(q(2,i))-L(4)*sin(q(1,i)).*cos(q(2,i)+q(4,i)) -L1_*cos(q(1,i)).*cos(q(2,i))-L(4)*cos(q(1,i)).*sin(q(2,i)+q(4,i)) -L(4)*cos(q(1,i)).*sin(q(2,i)+q(4,i));
%                 -L1_*cos(q(1,i)).*sin(q(2,i))+L(4)*cos(q(1,i)).*cos(q(2,i)+q(4,i)) -L1_*sin(q(1,i)).*cos(q(2,i))-L(4)*sin(q(1,i)).*sin(q(2,i)+q(4,i)) -L(4)*sin(q(1,i)).*sin(q(2,i)+q(4,i));
%                 0 -L1_*sin(q(2,i))+L(4)*cos(q(2,i)+q(4,i)) L(4)*cos(q(2,i)+q(4,i))];
%     invJ(:,:,i) = inv(invJ(:,:,i));
% end
% % size(invJ(:,:,1))
