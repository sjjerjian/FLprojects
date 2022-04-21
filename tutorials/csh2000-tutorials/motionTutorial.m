%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  File: motion-tutorial.m
%%%  Eero Simoncelli, 6/96.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This tutorial will take you through a simple derivation of an algorithm for
% computing optical flow, the projection of motion onto the image plane.  The
% algorithm is based on measurements of the gradient, but, as we will show, may
% also be thought of as a spatio-temporal "Energy" algorithm (ala
% Adelson/Bergen).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear
% The first portion of the tutorial will use only 1D examples.  These are
% easiar to visualize and easier to understand, and many of the important
% issues in the 2D case are also present in 1D.

% Consider a spot of light (an impulse), translating to the right at one pixel 
% per frame.  Motion is orientation.  We can view the sequence of signals in 
% the X (space) and T (time) dimensions, with amplitude corresponding to 
% intensity (time axis points down).

xtImpulse = upConv(eye(32),[1 1 1],'zero');
showIm(xtImpulse,'auto','auto','Intensity over space-time');

%%
% Now to compute motion, we can use the gradient method, which is
% based on the following constraint:
%
%      Ix . v + It = 0

% We'll need to construct derivative kernels separably from 1D filters:
%%gFilt = [0.233 0.534 0.233];
%%dFilt = [-0.459 0.0 0.459];
%gFilt = [0.0363 0.2489 0.4296 0.2489 0.0363];
%dFilt = [-0.1081 -0.2810 0.00 0.2810 0.1081];
X=[-5:5];  sig=2;
gFilt = exp(- (X.^2)/(2*sig^2));
dFilt = X .* gFilt;
subplot(2,1,1); lplot(gFilt);
subplot(2,1,2); lplot(dFilt);

%%
% Derivative filters in the X and T directions may be computed via separable
% combinations of the two 1D filters:

dxFilt = gFilt' * dFilt;
dtFilt = dFilt' * gFilt;
clf; showIm(dxFilt + i*dtFilt,'auto','auto',...
       'Space-time weighting functions of x- and t-deriv filters');

%%
% Derivatives in two dimensions have a beautiful rotation-invariance
% property: the derivative in an arbitrary directions is just a linear
% combination of the derivatives along the two axes:

angle = pi/6;

angleFilt = cos(angle) * dxFilt + sin(angle) * dtFilt;
showIm(angleFilt);

%%
% Because of this property, one  way to think about computing velocity
% is to find the angle at which we get a maximal filter response to
% the filter.  That is, first compute the filter responses in the
% center of the xt image:

filtSz2 = (size(dxFilt)-1)/2;
xtImpulseCtr = xtImpulse(16+[-filtSz2(1):filtSz2(1)], 15+[-filtSz2(2):filtSz2(2)]);
xFiltResp = sum(sum( dxFilt .* xtImpulseCtr));
tFiltResp = sum(sum( dtFilt .* xtImpulseCtr));

% ...and use the rotation-invariance property to compute the responses
% at a bunch of angles:
angles = (pi) * [-8:8]/8;
speeds = atan(angles);
responses = cos(angles) * xFiltResp + sin(angles) * tFiltResp;
clf; plot(angles, responses);
xlabel('filter preferred speed');
ylabel('response');


%%
% Instead of searching for the filter that gives the largest response,
% we can compute an estimate of velocity at each point in space and
% time using the gradient constraint:
dxImpulse = corrDn(xtImpulse,dxFilt,'dont-compute');
dtImpulse = corrDn(xtImpulse,dtFilt,'dont-compute');
vImpulse = - dtImpulse ./ dxImpulse;

% You should have gotten a "divide by zero" warning, because the dxImpulse
% image has lots of zeros.  Let's fix this with a hack (we'll fix it more
% properly in a moment):
dxImpulse = dxImpulse + (abs(dxImpulse)<1e-3);
vImpulse = - dtImpulse ./ dxImpulse;
clf; showIm(vImpulse,'auto','auto','Speed estimates over space-time');

% What you see is an x-t image of VELOCITY ESTIMATES.  Each point in the image
% corresponds to a different spatial position at a different time.  The
% intensity at each point is proportional to estimated speed.  look at a few
% values:
vImpulse(10:15,10:15)

%%
% The white pixels have a value of one, the correct speed.  But, there is an
% oblique line of zeroes in between the two lines of ones, and the remainder of
% the image is filled with zeros.

% The problem is that the spatial derivative dx-impulse is zero at some
% locations and we cannot compute the velocity at those points.  To "fix" the
% strip of zeroes between the lines of ones, we can combine the derivative
% measurements over a local neighborhood (by blurring the squared filter
% outputs).  THis comes from combining the derivative constraint over small
% neighborhoods, in a squared-error fashion:
%
%     E(v) = sum ( dx*v + dt )^2.
%
% The minimal solution for v is:
%
%     v = -sum(dx*dt) / sum(dx*dx)
%
% Note that we also add a small number (1e-4) to the demominator so that the
% quotient doesn't blow up.

% We'll need a blur kernel:
blurFilt = namedFilter('binom5')*namedFilter('binom5')';

dxImpulse = corrDn(xtImpulse,dxFilt,'dont-compute');
dtImpulse = corrDn(xtImpulse,dtFilt,'dont-compute');

vBlurImpulse = - (corrDn(dtImpulse*dxImpulse,blurFilt,'dont-compute')) ./ ...
                   (1e-4 + (corrDn(dxImpulse*dxImpulse,blurFilt,'dont-compute')));
vBlurImpulse(1:6,1:6)
clf; showIm(vBlurImpulse,'auto','auto','Speed estimates over space-time');

%%
% The velocity estimate is now correct (v=1) in the vicinity of the moving
% impulse.  You still get zero far from the diagonal This is because there is NO
% motion information in the portion of the images away from the impulse.  The
% computation is thus singular at these points.  As we will see, this problem
% occurs in an even more serious way in 2D.

% To save typing, we've packaged this mini-algorithm into a function called
% "compute1dFlow".

clf; showIm(compute1dFlow(xtImpulse),'auto','auto',...
       'Speed estimates over space-time');

% Now we'll compute velocity on a one-dimensional image with more interesting
% content.  We'll make a 1D fractal noise pattern moving one pixel to the left
% per time-step:
fract = mkFract([1 32],2);
clf; plot(fract);
xlabel('Position')
ylabel('Intensity (arbitrary units)')

xtFract  = zeros(16,size(fract,2));
for row = 1:size(xtFract,1)
	xtFract(row,:) = [fract(1,(33-row):32), fract(1,1:(32-row))];
end

vBlurFract = compute1dFlow(xtFract);
clf; showIm(vBlurFract,'auto','auto','Speed estimates over space-time');
drawnow

% Now the entire space-time image of speed estimates contains values close
% to 1 because the entire pattern is moving 1 pixel per frame.

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Motion in the frequency domain: a connection between the gradient and
% spatio-temporal energy mechanisms.

% A translating one-dimensional pattern has a fourier spectrum lying on a line
% through the origin.  The slope of the line corresponds to the pattern
% velocity.

fractMag = fftshift(abs(fft2(xtFract)));
clf; showIm(fractMag,[0 100],'auto','Spatiotemporal Fourier amplitude');
drawnow

% The Fourier magnitude image is plotted over a range of [-pi, pi] in both the
% wx (spatial frequency) and wt (temporal frequency) directions.  Most of the 
% energy in this stimulus lies along an oblique line in wx-wt.  
% [But why doesn't it lie EXACTLY on a line?]

%%
% Remember that convolution corresponds to multiplication in the frequency
% domain.  Therefore, the Fourier Transform of the derivative images (dxFract
% and dtFract) corresponds to the product of the DFT's of the filter and
% fract.  Let's look at the Fourier magnitude of the filters (first the time
% derivative filter, the x-derivative filter):

dtMag = fftshift(abs(fft2(dtFilt,32,32)));
dxMag = fftshift(abs(fft2(dxFilt,32,32)));
clf; showIm(dxMag + i*dtMag,'auto','auto',...
       'Spatiotemporal freq responses of t- and x-deriv filters');
drawnow

%%
% The derivative operators are spatiotemporal filters.  Not only that, we can
% write the temporal derivative as a sum of two other spatio-temporal filters
% oriented at +/- 45 degrees (that is, filters that are most sensitive to
% either rightward or leftward motion):

drFilt = (dxFilt+dtFilt)/2;
dlFilt = (dxFilt-dtFilt)/2;
clf; showIm(dlFilt + i*drFilt,'auto','auto',...
       'Space-time weighting functions of right- and left-selective filters');
drawnow

%%
% Now look at the FTs of these two direction selective filters.  They are
% rightward and leftward spatiotemporal energy filters.

rMag = fftshift(abs(fft2(drFilt,32,32)));
lMag = fftshift(abs(fft2(dlFilt,32,32)));
clf; showIm(rMag + i*lMag,'auto','auto',...
       'Spatiotemporal freq responses of right and left-selective filters');
drawnow

%%
% Given this relationship between right-/left-selective filters and the
% x-/t-derivative filters, we can now rewrite the velocity calculation to 
% look like a standard Adelson-Bergen style spatio-temporal energy calculation, 
% ((R - L)/S).

dlFract = corrDn(xtFract,dlFilt,'dont-compute');
drFract = corrDn(xtFract,drFilt,'dont-compute');
dxFract = corrDn(xtFract,dxFilt,'dont-compute');

lEnergy = corrDn(dlFract.*dlFract, blurFilt,'dont-compute');
rEnergy = corrDn(drFract.*drFract, blurFilt,'dont-compute');
sEnergy = corrDn(dxFract.*dxFract, blurFilt,'dont-compute');
sEnergy = sEnergy + 1e-4;  %% avoid divide-by-zero

vBlurSteFract = (lEnergy - rEnergy) ./ sEnergy;

% Print statistics to show that these two things are identical (up to 1e-15
% or so).
imStats(vBlurFract,vBlurSteFract);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Temporal Aliasing:

% The patterns above were moving only one pixel per frame.  Now consider a
% pattern moving two pixels per frame:
fract2 = mkFract([1 64],2);
xtFract2  = zeros(16,size(fract2,2));
for row = 1:size(xtFract2,1)
	xtFract2(row,:) = [fract2(1,(65-row*4):64), fract2(1,1:(64-row*4))];
end
clf; showIm(xtFract2,'auto','auto','Intensity over space-time');
drawnow

%%
% Now let's try to estimate the velocity of this stimulus:
vBlurFract2 = compute1dFlow(xtFract2);
clf; showIm(vBlurFract2,'auto','auto','Speed estimates over space-time');
drawnow

%%
% Look at a histogram of the estimated velocities:
clf; hist(vBlurFract2(:),30);

%%
% The speed estimates have quite a lot of variability.  What is happening?  The
% filters are too small for an image displacement of 2 pixels, so the algorithm
% starts to fail.  To see this another way, look at the Fourier transform of
% the signal:
fract2Mag = fftshift(abs(fft2(xtFract2,32,32)));
clf; showIm(fract2Mag,[0,100],'auto','Spatiotemporal Fourier amplitude');
drawnow

% You should be able to see aliased copies of the spectrum, parallel to the
% central "stripe".  Why does this happen?  How does it affect the derivative
% filters (look back at dxMag and dtMag)?

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% A Coarse-to-fine algorithm:

% One way to get around the temporal aliasing problem is to use filters of
% lower spatial frequency response.  These will not "see" the aliased copies,
% and thus will be able to get a better velocity estimate.  Alternatively, we
% can blur and subsample the images spatially, thus reducing the motion
% displacement per frame.

xtFract2Blur = corrDn(xtFract2,namedFilter('gauss5')','circular',[1 2]);
clf; showIm(xtFract2Blur,'auto','auto','Coarse scale intensity over space-time');
drawnow

%%
% Note that the aliasing in the power spectrum will now have a much smaller
% effect on the filters because most of the energy is now back on the main
% diagonal:

xtFract2BlurMag = fftshift(abs(fft2(xtFract2Blur,32,32)));
clf; showIm(xtFract2BlurMag,[0,100],'auto','Spatiotemporal Fourier amplitude');
drawnow

%%
% Now we can compute motion on this blurred and subsampled image.
% Rember that the resulting velocities must now be multiplied by 2
% to get them in units of pixels/frame relative to the sampling of
% the original images.

v2BlurFract2 = 2 * compute1dFlow(xtFract2Blur);
clf; showIm(v2BlurFract2,'auto','auto','Speed estimates');
drawnow

% We could stop here. But the current velocity field is computed at reduced
% resolution, and we may want higher resolution velocity fields.  In
% particular, for real images, the motion is often NOT uniform translation.
% Thus we want to use the smallest filters possible to compute the most LOCAL
% velocity possible.

% One solution is to use a "coarse-to-fine" procedure.  We get an initial
% estimate at the coarse scale using big filters (or equivalently, using
% blurred and subsampled images).  Then we refine that estimate at successively
% finer scales.  After computing the initial coarse scale estimate of the
% velocity, we try to "undo" motion by aligning the images according to the
% coarse-scale estimate.  Then the finer scale filters are used to compute a
% correction to the coarse-scale flow field.

% To undo the motion, we translate or "warp" each frame (each scan line of the
% X-T diagram) back toward the center frame according to our current motion
% estimate.

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 2 Dimensional Images.

% The velocity estimation problem becomes a bit more difficult, both
% conceptually and computationally, when considering 2D images.  Nevertheless,
% the same basic gradient constraint may be used as a basis for a reasonable
% algorithm....

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
