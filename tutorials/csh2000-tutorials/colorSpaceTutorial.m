% colorSpaceTutorial
%
% This tutorial introduces some ideas about
% color space transformations.
%
% This tutorial was created by updating an existing
% tutorial written by Brian Wandell.  It builds on
% ideas developed in the colorRenderingTutorial, so
% you should go through and understand that first.
%
% See also: colorRenderingTutorial
%
% Dependencies: a) imshow (image processing toolbox)
%               b) Psychophysics Toolbox:PsychColorimetricData
%               b) colorFunctions subfolder
%
% 06/18/98  dhb   Updated Wandell version for CSH '98.
%				  Removed stuff redundant with colorRenderingTutorial.
%                 Made data structures to match Psychophysics Toolbox
%                 conventions.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Introduction
%
% In this tutorial, we will introduce some basic color space
% transformations.
%
% You will notice that much of the work below involves matrix
% algebra.  We suggest that you pull out a pad of paper and draw
% matrix tableaus as we proceed to help you follow the calculations.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Initialize.
%
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Transformation between two color spaces, each specified
% by its own set of color matching functions.  Exact case.
%
% Often we want to transform data between two color spaces.
% For example we might have CIE XYZ tristimulus coordinates
% and want to convert to cone coordinates.
%
% The easiest case to handle is when the two color spaces
% represent the same standard observer, that is, the two sets
% of color matching functions are related by an 3 by 3 linear
% transformations.  For example, this is the case for the
% Judd-Vos XYZ functions and the Smith-Pokorny estimates
% of the cone fundaments.  Lets verify that.

% Load the data.  If you want to know more about the format
% in which the data are stored, type "help PsychColorimetricData"
load T_cones_sp
load T_xyzJuddVos

% In a better world, we would have called the file that
% contained the Judd-Vos XYZ functions T_XYZJuddVos, and
% similarly for the file that contains the 1931 XYZ functions
% which are loaded later in this tutorial.  This would avoid
% confusion between tristimulus (XYZ) values and chromaticity
% (xyz) values.  However, the naming convention is used
% widely in other code and we decided not to change it for
% this tutorial.

% Figure out spectral sampling.  Vector S_cones_sp specifices
% wavelength sampling in form [start delta number] and this
% code produces actual wavelength values.
S = S_cones_sp;
start = S(1); delta = S(2); number = S(3);
spectrum = (start:delta:(start+(number-1)*delta))';

% Plot each set of functions
figure(1); clf;
plot(spectrum,T_cones_sp(1,:)','r');
hold on
plot(spectrum,T_cones_sp(2,:)','g');
plot(spectrum,T_cones_sp(3,:)','b');
xlabel('Wavelength (nm)');
ylabel('Relative Sensitivity');
title('Cone Spectral Sensitivity Functions');
hold off
figure(2); clf;
plot(spectrum,T_xyzJuddVos(1,:)','r');
hold on
plot(spectrum,T_xyzJuddVos(2,:)','g');
plot(spectrum,T_xyzJuddVos(3,:)','b');
xlabel('Wavelength (nm)');
ylabel('Tristimulus value');
title('Judd-Vos Color Matching Functions');
hold off
drawnow;

% Use regression to find best 3 by 3 linear transformation.
% See MATLAB's "help mldivide" for information on what the magic
% \ operator does.
M_JuddVosToCones = ((T_xyzJuddVos')\(T_cones_sp'))';

% Use the matrix to compute cone sensitivity from Judd-Vos.
% There should be essentially perfect agreement.  To see
% this, we plot the estimates over the original curves.
T_cones_fromJuddVos = M_JuddVosToCones*T_xyzJuddVos;
figure(1);
hold on
plot(spectrum,T_cones_fromJuddVos(1,:)','r*');
plot(spectrum,T_cones_fromJuddVos(2,:)','g*');
plot(spectrum,T_cones_fromJuddVos(3,:)','b*');
hold off
drawnow;

% We can now apply the transformation matrix to any
% Judd-Vos tristimulus coordinates to get cone
% coordinates.  To have something to apply to,
% we'll compute the Macbeth color checker coordinates
% with the Judd-Vos color matching functions and then
% transform them to cone space.  The computation
% is done out in detail in the colorRenderingTutorial.
load sur_macbeth ; load spd_D65
spectralSignals = diag(spd_D65)*sur_macbeth;
xyzJuddVos = T_xyzJuddVos*spectralSignals;
coneSignals = T_cones_sp*spectralSignals;
coneSignals_fromJuddVos = M_JuddVosToCones*xyzJuddVos;

% Make a little plot to verify that the transformation
% worked.  Each panel plots (for one cone coordinate)
% the value obtained by transforming the XYZ values vs
% the dirctly computed cone value.  If the transformation
% is correct, all points should lie on the positive diagonal.
figure(1); clf;
subplot(1,3,1);
plot(coneSignals(1,:),coneSignals_fromJuddVos(1,:),'r+');
plotMax = max(max([coneSignals(1,:),coneSignals_fromJuddVos(1,:)]));
axis([0 plotMax 0 plotMax]);
axis('square');
hold on
plot([0 plotMax]',[0 plotMax]','r');
hold off
xlabel('Direct L');
ylabel('Transform');
subplot(1,3,2);
plot(coneSignals(2,:),coneSignals_fromJuddVos(2,:),'g+');
plotMax = max(max([coneSignals(2,:),coneSignals_fromJuddVos(2,:)]));
axis([0 plotMax 0 plotMax]);
axis('square');
hold on
plot([0 plotMax]',[0 plotMax]','g');
hold off
xlabel('Direct M');
subplot(1,3,3);
plot(coneSignals(3,:),coneSignals_fromJuddVos(3,:),'b+');
plotMax = max(max([coneSignals(3,:),coneSignals_fromJuddVos(3,:)]));
axis([0 plotMax 0 plotMax]);
axis('square');
hold on
plot([0 plotMax]',[0 plotMax]','b');
hold off
xlabel('Direct S');
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Transformation between two color spaces, each specified
% by its own set of color matching functions, but where the
% two are not exact linear transformations.
%
% We can apply the same method, but the result will not
% be exact.  The code below repeats the code above, but we
% will use the 1931 CIE XYZ color matching functions, rather
% than the Judd-Vos.  The 1931 XYZ color matching functions
% are not an exact linear transformation away from the Smith
% Pokorny fundamentals, so the two represent different observers.

% Calculate transformation, just as before.
load T_xyz1931
M_1931ToCones = ((T_xyz1931')\(T_cones_sp'))';

% Compare the cone sensitivites derived from 1931 XYZ
% to Smith-Pokorny estimates.  Note that now the match
% is not exact.  This is because the two color spaces
% represent different standard observers (i.e. they
% are not within a linear transformation of one another).
% The mismatch is mainly for the S-cones, and this is
% because it is in the short wavelengths that the
% Judd-Vos corrections were applied.
T_cones_from1931 = M_1931ToCones*T_xyz1931;
figure(1); clf
plot(spectrum,T_cones_sp(1,:)','r');
hold on
plot(spectrum,T_cones_sp(2,:)','g');
plot(spectrum,T_cones_sp(3,:)','b');
xlabel('Wavelength (nm)');
ylabel('Relative Sensitivity');
title('Cone Spectral Sensitivity Functions');
plot(spectrum,T_cones_from1931(1,:)','r*');
plot(spectrum,T_cones_from1931(2,:)','g*');
plot(spectrum,T_cones_from1931(3,:)','b*');
hold off
drawnow;

% We can see the same fact by comparing coordinates
% for the Macbeth color checker.  This parallels
% what was done above.
xyz1931 = T_xyz1931*spectralSignals;
coneSignals_from1931 = M_1931ToCones*xyzJuddVos;
figure(1); clf;
subplot(1,3,1);
plot(coneSignals(1,:),coneSignals_from1931(1,:),'r+');
plotMax = max(max([coneSignals(1,:),coneSignals_from1931(1,:)]));
axis([0 plotMax 0 plotMax]);
axis('square');
hold on
plot([0 plotMax]',[0 plotMax]','r');
hold off
xlabel('Direct L');
ylabel('Transform');
subplot(1,3,2);
plot(coneSignals(2,:),coneSignals_from1931(2,:),'g+');
plotMax = max(max([coneSignals(2,:),coneSignals_from1931(2,:)]));
axis([0 plotMax 0 plotMax]);
axis('square');
hold on
plot([0 plotMax]',[0 plotMax]','g');
hold off
xlabel('Direct M');
subplot(1,3,3);
plot(coneSignals(3,:),coneSignals_from1931(3,:),'b+');
plotMax = max(max([coneSignals(3,:),coneSignals_from1931(3,:)]));
axis([0 plotMax 0 plotMax]);
axis('square');
hold on
plot([0 plotMax]',[0 plotMax]','b');
hold off
xlabel('Direct S');
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Chromaticity coordinates.
%
% The standard CIE xy chromaticity diagram is a two
% dimensional representation of the XYZ tristimulus
% values.  In going from three dimension to two,
% information is lost.  In particular it is not
% possible to recover the overall intensity of a
% stimulus, nor its luminance, from its chromaticity
% representation.
%
% The formulae for computing chromaticity coordinates
% are simple:
%   x = X/(X+Y+Z);
%   y = Y/(X+Y+Z);
%
% Sometimes people define a little z chromaticity as
%   z = Z/(X+Y+Z);
% but this adds no information since z = 1 - x - y.
%
% As an example, we will compute the chromaticity coordinate
% for the monochromatic lights specified by the variable
% spectrum.  Note that each column of T_xyz1931 (3 by 81)
% contains the XYZ coordinates of a monochromatic light.
spectral_x = T_xyz1931(1,:) ./ sum(T_xyz1931);
spectral_y = T_xyz1931(2,:) ./ sum(T_xyz1931);

% Make a plot
figure(1); clf;
plot(spectral_x',spectral_y','k');
axis('square');
axis([0 1 0 1]);
drawnow;

% For fun, let's now add the chromaticities of a monitor's
% phosphors to our plot.
load B_monitor
monitorXYZ = T_xyz1931*B_monitor;
monitor_x = monitorXYZ(1,:) ./ sum(monitorXYZ);
monitor_y = monitorXYZ(2,:) ./ sum(monitorXYZ);
figure(1);
hold on
plot([monitor_x monitor_x(1)],[monitor_y monitor_y(1)],'k');
plot(monitor_x(1),monitor_y(1),'r*');
plot(monitor_x(2),monitor_y(2),'g*');
plot(monitor_x(3),monitor_y(3),'b*');

% The triangle in the figure shows the gamut of
% chromaticites that can be produced on this
% monitor.  Can you prove this?


