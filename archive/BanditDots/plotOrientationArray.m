% plotOrientationArray

clear

% mean and stdev of the distribution of orientations in the array
mu = 0; % 0 is vertical
sigma = 16; 
nLines = 50; % number of line segments in the array
len = 0.1; % length of each segment

ori = normrnd(mu,sigma,nLines,1); % draw random orientations


figure(1); clf;
for n = 1:nLines
    
    r = len;
    theta = ori(n)+90; % rotate 90 deg so that 0 is vertical
    
    X = r*cosd(theta);
    Y = r*sind(theta);
    
    % randomize origin of each line
    origX = rand; origY = rand;
    
    plot([origX origX+X],[origY origY+Y],'k-'); hold on;
    
end