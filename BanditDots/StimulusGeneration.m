% plotOrientationArray

clear

% mean and stdev of the distribution of orientations in the array
mu = 0; % 0 is vertical
sigma = 16; % sigma is jitter
nLines = 50; % number of line segments in the array
len = 0.1; % length of each segment
% max x/y coordinate is always 1 + len

circlecenter = [(1+len)/2, (1+len)/2];
ccx = (1+len)/2;
ccy = (1+len)/2;
circleradius = (1+len)/4;

ori = normrnd(mu,sigma,nLines*10,1); % draw random orientations
countInside = 0;
n=1;

figure(1); clf;
while countInside < nLines && n< nLines*7
    n= n+1;
    r = len;
    theta = ori(n)+90; % rotate 90 deg so that 0 is vertical
    
    X = r*cosd(theta);
    Y = r*sind(theta);
    
    % randomize origin of each line
    x1 = rand; y1 = rand;
    x2= x1+X;
    y2 = y1+Y;
    if inside(x1,y1,x2,y2,circleradius,ccx,ccy)==true
        plot([x1 x2],[y1 y2],'k-'); hold on;
        countInside=countInside+1;
    end
end

axis([0 1.2 0 1.2])
axis equal;
disp("Lines inside circle are: " + countInside);
viscircles(circlecenter,circleradius,'color','none')


function in = inside(x1,y1,x2,y2,circleradius,ccx,ccy)
    if ((x1-ccx)^2 + (y1-ccy)^2) > circleradius^2
        in = false;
    elseif ((x2-ccx)^2 + (y2-ccy)^2) > circleradius^2
        in = false;
    else
        in = true;
    end
end