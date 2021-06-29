function yq = linear_pcw_fit(y,x,step)

% linear piece-wise fit
% effectively, for vals with x-coordinate x, interpolate between
% each one, with stepsize step

% x(1) and x(end) are the ends of the data

xq = x(1):step:x(end);
yq = nan(size(xq));

start = 1;
for ix=1:length(x)-1
    
    [~,stop] = min(abs(xq-x(ix+1)));
    
    yq(start:stop) = linspace(y(ix),y(ix+1),numel(xq(start:stop)));
   
    start = stop+1;
end
    

