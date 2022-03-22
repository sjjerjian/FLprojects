%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulating a Differential Equation.

% We want to solve (numerically) the following differential
% equation:
%
%    tau dy/dt + y(t) - x(t) = 0
%
% where x(t) is the input, y(t) is the output, and tau is the
% time constant that determines the behavior of the system
%
% To do this, we discretize time and rewrite the equation as a
% finite difference equation:
%
%    tau deltaY[i]/deltaT + y[i] - x[i] = 0
%
% where i is an integer index that corresponds to discrete time
% intervals, and where
%
%    y[i+1] = y[i] + deltaY[i]
%    
% Then with a little bit of algebra, we write:
%
%    y[i+1] = y[i] + (deltaT/tau) (x[i] - y[i])

% First we need to define deltaT (the time step)
deltaT=1e-3;				% 1e-3 secs = 1 msec

% Next we need to define the range of times for the simulation,
% and the indices for the discrete time samples:
times=[0:deltaT:1];			% 1 sec simulation
indices=[1:length(times)];

% Next we need to define a discretely sampled input.  I'll make
% two different inputs, one a sinusoid, and the other a step.
% First, the step.

% Initialize it to zeros:
step=zeros(size(indices));
% Set values after the first 1/4 sec to ones:
step(times>.25)=ones(size(step(times>.25)));
% Plot it:
plot(times,step);
xlabel('Time')
ylabel('Input')

% Next the sinusoid:
frequency=4; 				% 4 Hz
sine=sin(2*pi*frequency*times);
plot(times,sine);
xlabel('Time')
ylabel('Input')

% Now that we have some inputs to work with, we need to specify
% the time constant parameter.
tau=0.1					% 100 msec

% Finally, we need to specify the initial condition:
y(1)=0;

% Then we choose one of the two inputs, and iterate through time,
% updating y as we go:
x=step;
for i=[1:length(times)-1]
  y(i+1) = y(i) + (deltaT/tau) * (x(i) - y(i));
end
% Plot it:
plot(times,y);
xlabel('Time')
ylabel('Output')

% You can plot a segment of the output with:
segment=[240:260];
plot(times(segment),y(segment));
xlabel('Time')
ylabel('Output')

% Here are some things for you to do:
%
% - Try the sinusoidal input.
%
% - Try changing the time constant tau to make it a bit longer or
%   a bit shorter.
%
% - Try making tau really short?  Specifically, set tau to .5
%   msec and run the step input.  Then try making tau even shorter,
%   e.g., tau=.1 msec.  Why does the simulation behave this way?
%
% - Now fix the time constant at some reasonable value, and run a
%   bunch of sinusoidal inputs with different frequencies.  What
%   happens and why?

