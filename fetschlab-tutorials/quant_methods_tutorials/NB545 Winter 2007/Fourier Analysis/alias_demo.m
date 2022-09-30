% Aliasing demo
% Written for class, Jan 23, 2003
% M. N. Shadlen
% This is minimally documented right now, but it could be converted into a
% tutorial

% Define a discrete time axis. The spacing is dt
dt = .0001;      % our units are seconds. So each integer time point is scaled to dt sec
tmax = 1
t = [0:dt:tmax-dt]';
% Notice that we are sampling time at a rate of 1/dt. This will be
% important later.
samplingRate = 1./dt    % in Hz

% The original sample I'm going to play with is an amplitude modulated (AM)
% sinusoid. I use this to entice some students to ask questions about
% missing fundamental illusion and some radio basics. It could promote
% confusion too. We'll see.
s1 = 0.5*(1 + cos(2*pi*110*t)) .* cos(2*pi*880*t);
% Another function that is fun to see or think about
% s1 = 0.5*(1 + cos(2*pi*110*t)) .* rand(size(t));

figure(1), clf
plot(t,s1)
newSamplingRate = 100;
comb = mod(t,1./newSamplingRate)==0;
sum(comb)       % check to make sure the modular arithmetic worked
% Create undersampled signal
s1u = comb .* s1;

% Listen to the waveforms
sound(s1,samplingRate)
sound(s1u,samplingRate)

nyq = samplingRate/2;  % nyquist rate
fax = [0:nyq]'; % create frequency axis (positive frequencies only)
n = length(fax);
% Show time domain on the left, frequency domain on right
figure(4), clf
clear ax
ax(1) = subplot(3,2,1);
plot(t,s1)
ax(2) = subplot(3,2,2);
s1 = abs(fft(s1));
plot(fax,s1(1:n))
ax(3) = subplot(3,2,3);
plot(t,comb)
ax(4) = subplot(3,2,4);
COMB = abs(fft(comb));
plot(fax,COMB(1:n))
ax(5) = subplot(3,2,5);
plot(t,s1u)
ax(6) = subplot(3,2,6);
S1U = fftshift(abs(fft(s1u)));
plot(fax,S1U(1:n))
set(ax(1:2:6),'XLim',[0 .2])
set(ax(2:2:6),'XScale','Linear','XLim',[1 2000])

% Look at frequency on log scale
set(ax(2:2:6),'XScale','Log','XLim',[1 2000])


% Another example. THIS IS BROKEN. I'll fix it one of these days.

% Let's illustrate with a more realistic set of functions.
s1 = sin(2*pi*3*t) + .33 * sin(2*pi*9*t) + .2 * sin(2.*pi*15*t);
% Corrupt this signal by adding noise
s1 = s1 + rand(size(s1));
% You could use a different blurring function. I'll apply a gaussian blur.
% See if you can figure out what the fftshift is doing here.
g = fftshift(exp(-((t-.5)/t).^2));
s2 = fft(g) * fft(s1);
s2 = ifft(s2);
figure(4), clf
ax(1) = subplot(3,2,1)
plot(t,s1)
ax(2) = subplot(3,2,2)
s1 = abs(fft(s1));
plot(fax, s1(1:n));
ax(3) = subplot(3,2,3)
plot(t,g)
ax(4) = subplot(3,2,4)
G = abs(fft(g));
plot(fax,G(1:n))
ax(5) = subplot(3,2,5)
plot(t,s2)
ax(6) = subplot(3,2,6)
plot(fax, abs(s2(1:n)));
set(ax(1:2:6),'XLim',[0 .2])
set(ax(2:2:6),'XScale','Linear','XLim',[1 2000])




