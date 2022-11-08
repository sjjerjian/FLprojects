% extraFourierStuff

dt = .01;      % our units are seconds. So each integer time point is scaled to dt sec
tmax = 1
t = [0:dt:tmax-dt]';
% Notice that we are sampling time at a rate of 1/dt. This will be
% important later.
samplingRate = 1./dt    % in Hz

nyq = samplingRate/2

k = [0:nyq];
% Make a matrix representing the cosine basis set
[T,K] = ndgrid(t,k);

C = cos(2*pi*K.*T);
S = sin(2*pi*K.*T);
imagesc(C(:,[1:10]))
colormap(gray)
set(gca,'Visible','off')
print -dtiff C1thru10.tif
imagesc(S(:,[1:10]))
colormap(gray)
set(gca,'Visible','off')
print -dtiff S1thru10.tif

clf
imagesc(C)
set(gca,'Visible','off')
print -dtiff C_all.tif







sqwav = sign(cos(3*pi*2*t)+eps);
q = fft(sqwav);
a = real(q(1:nyq+1));
b = imag(q(1:nyq+1));
max(abs(a))
max(abs(b))

clf
imagesc(sqwav)
set(gca,'Visible','off')
pos = get(gca,'Position')
set(gca,'pos',pos .* [1 1 .1 1])
print -dtiff sq3.tif

clf
imagesc(a)
set(gca,'Visible','off')
set(gca,'pos',pos .* [1 1 .1 1])
print -dtiff fftReal_sq3.tif


B = repmat(a',length(t),1) .* C;
clf 
imagesc(B)
set(gca,'Visible','off')
print -dtiff cosDecomp_sq3.tif

sum(a>1)
plot(a)
a(4)
alist = [4:6:length(a)]
n = length(alist)
psum = zeros(size(B(:,1)));
clf
i=1
for i = 1:n
    f = B(:,alist(i));
    ax(2*i - 1) = subplot(n,2,2*i-1);
    plot(t,f)
    ax(2*i) = subplot(n,2,2*i);
    psum = psum+f;
    plot(t,psum);
end
set(ax,'Box','off','YLim',[-100 100])
set(ax,'YTick',[])
set(ax(1:end-2),'XTickLabel',[])
set(ax,'TickDir','out')



    



