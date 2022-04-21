% leastsqTutorial.m  script
%
% 
% This is a beginner's tutorial.  If you play with the concepts in this 
% tutorial, you can go from knowing next to nothing about linear algebra to an
% intuition that will help with multiple linear regression. The approach is
% based on Strang's 'Linear Algebra'. You will learn about subspaces
%
% M.N. Shadlen for Cold Spring Harbor 1998
% Aug 2004, improvements in nomenclature
% Mar 2005, corrected typos in lines 71-73 (Thanks Adrienne!)
% Feb 2009, organizing by cells. Minor clean up

%% The model problem
clear all
% suppose you have ordered pairs
data = [-1 1 2; 1 1 3]'

% Our job is to find the best fitting line for these three points.
% Here are the points
hf1 = figure(1), clf
plot(data(:,1), data(:,2),'o')
set(gca,'XLim',[-3 3],'YLim',[-1 4])

%% The matrix approach
% We are going to set this up as a matrix equation.  y = A*b. A contains our
% observations, which in this are the x values, along with a column of
% ones. The matrix, A, is called the design matrix for reasons that will
% be evident later. Y represents the dependent variable, what we
% typically graph on the y-axis. In one sense these are observations or
% measurments. In another sense, they are predictions. I will try to be
% clear about this distinction: y_obs is a vector of observed values for
% y. Depending on context, y or y_pred or y_guess is the predicted value. The
% remaining variable, b, which is also a vector, contains the fit
% coefficients. In the case of a line, there are just two: intercept and
% slope.

% [Aside for readers of Strang. Strang uses x where I'm using b. That's
% reasonable, because this is the unknown vector that we need to solve for.
% Most scientists naturally think of x as the argument of a function and
% they customarily represent a line as y = b0 + b1*x. To make things easier,
% I'm using this more common notation. In a moment, we'll see that x (in my
% notation, not Strang's) is just a column of A.]

A = [ones(3,1) data(:,1)]

% Notice that the first column of A is just a column of ones. 
% The second column contains the values plotted on the abscissa -- what
% we would normally call the x-variable: the values -1 1 2.

x = data(:,1)

% The observed y values are the 2nd column of the data
y_obs = data(:,2)

% We know A and we know x.  We want to find the best solution for b. 
% The first thing you need to convince yourself of is that the 2 by 1 column 
% vector, b, contains the constant and slope of the best fitting line.
% This is a good time to remind yourself of how a matrix multiplies a vector.

% Let's take a wild guess at b.  Let's guess that the intercept is .5 and
% the slope is 1
b = [.5 1]'

% Notice that you get out 3 values from A*x
y_guess = A*b
hold on, plot(x,y_guess,'+',x,y_guess,':')

% That's a pretty crummy fit, but that's not the point.  All you're supposed
% to see right now is that by going across the rows of A and down the
% column vector x, you get the equation 
[	A(1,1)*b(1) + A(1,2)*b(2)
	A(2,1)*b(1) + A(2,2)*b(2)
	A(3,1)*b(1) + A(3,2)*b(2)
]


% If you don't understand that last equation, stop.  There's no point
% in going any further. Ask for help or pick up a book.  
% If you understand the equation, then ask yourself the following.
% Suppose we made A as a column of 1s a column like A(:,2) and a third
% column of the squares:  A(:,3) = A(:,2).*A(:,2).  We would be looking for 
% a column vector, b, that has 3 components.  Do you see that these are the
% coefficients of the best fitting quadratic? b1 + b2*x + b3*x^2 = y
% Best to convince yourself of that too before going on.

% Now this way of setting up the equation jives with your sense of matrix
% times vector by going along the rows of the matrix and multiplying the
% elements by those going down the column vector.  Another way to put that
% is that each element in y is a dot-product of the corresponding row in A
% with the column vector b.  That's worth thinking about, but it does not
% lend itself to a particularly lucid view of the regression problem.  For
% this, it is worth thinking about the matrix times the vector in a
% different way.  Instead of going along the rows of A, consider the
% columns.  This is actually the easier way to think about it.  Notice that
% b is the same dimension as a column in A (3 by 1).  Convince yourself that
% y is just the weighted sum of the columns in A.  The weights come from the
% elements of b.

y_guess = b(1)*A(:,1) + b(2)*A(:,2)


% Again, it's not worth going on if you don't see this version of A*b

% Now here is where things start getting interesting.
% Once we recognize that y is the weighted sum of just 2 vectors, 
% it must be the case that y lies in the plane that these vectors 
% define.  No matter how you cut it, you can only find vectors, y, that you 
% can get to by adding the two vectors that form the columns of A.
% This plane is called the column-space of A.  It includes the origin.

%% 3D Picture
% It helps to look at a picture.  The way I'm going to draw this picture
% includes a step or two that we're not ready for.  Just execute it
% for the time being.  Later, we'll come back to why the matrix algebra
% works.


% Let's start by drawing our two column vectors
hf2 = figure(2)
hold off; hp = plot3(A(1,:),A(2,:),A(3,:),'ro','MarkerFaceColor','r'), hold on;
hl1= line([0 A(1,1)], [0 A(2,1)], [0 A(3,1)], 'Color', 'r')
hl2=line([0 A(1,2)], [0 A(2,2)], [0 A(3,2)],'Color', 'r')
s = sprintf(' text(%d,%d,%d,''(%d,%d,%d)'')', [A(:,1) A(:,1)]), eval(s)
s = sprintf(' text(%d,%d,%d,''(%d,%d,%d)'')', A(:,[2 2])), eval(s)
set(gca,'Box', 'on','XLim',[-4 4],'YLim',[-4 4],'ZLim',[-2 3])
set(gca,'XGrid','on','YGrid','on','ZGrid','off')

% You can imagine the plane
hfill = fill3([0 A(1,:) 0], [0 A(2,:) 0], [0 A(3,:) 0], 'b','EraseMode','xor')

% Now here's the part that you need to ignore for a moment.
% The next line makes a projection matrix. We're going to use it 
% to project the xy plane into the column space of A.
P = A * inv(A'*A) * A'
q = P * 2 * [1 1 1; 1 -2 1; -2 -2 1; -2 1 1; 1 1 1]'
h2 = fill3(q(1,:), q(2,:), q(3,:) , 'b','EraseMode','xor','EdgeColor','r')
delete(hfill)
set(gca,'FontSize',14)
xlabel('dim_1'), ylabel('dim_2'), zlabel('dim_3')

figure(2)
for k = 332:-4:300
    set(gca,'View',[k 30])
			drawnow
end
for k = 30:5:50
  set(gca,'View',[300 k])
	drawnow
end


% Are you convinced that no matter what we guess for our solution, x, we
% will end up with something in the plane?  You should have convinced yourself
% of this earlier.  But here is a quick demo.


% It is important to realize that any guess we make for the answer to
% our problem will make a vector in the plane spanned by A's column vectors.
% We'll make 5 random guesses for g and see what we get out
for i = 1:5
	g = 2 * (rand(2,1) - .5)
	res = A * g
	hg(i) = plot3(res(1),res(2),res(3),'c.')
	hl(i) = line([0 res(1)], [0 res(2)], [0 res(3)], 'Color','y')
end

% You might want to replay the axis-rotation again. Depending on where
% those random vectors point, it may or may not be obvious that they are in 
% the plane.

% Now that you have convinced yourself that no matter what we choose for x,
% we are going to get out (after multiplication by A) a vector that lies in
% the plane.  That is unfortunate because the perfect solution would be the
% vector y_obs, and it is usually the case that y_obs does not lie in our
% plane!  What a bummer.  But then again, if it did, we would have a perfect
% solution -- and that's not what we're dealing with in least squares.

% Let's look at the vector y_obs.  Remember this is the list of ordinate
% values for the points in Figure 1.

plot3(y_obs(1), y_obs(2), y_obs(3),'g*')
line([0 y_obs(1)], [0 y_obs(2)], [0 y_obs(3)], 'Color','g')
set(gca,'ZLim',[-1 3])

% It helps to get rid of the extra lines
delete(hl)
delete(hg)

% It also helps to anchor the orgin 

line([0 0],[0 0],[-1 0],'LineWidth',3,'Color','r')

%% a little spinning couldn't hurt
figure(2)
k = 50
i = 300
set(gca,'View',[i k])
drawnow, % pause
for k = 50:-4:30
 	set(gca,'View',[i k])
	[i k]
	drawnow
	% pause
end
for i = 300:-5:220
 	set(gca,'View',[i k])
	[i k]
	drawnow
	% pause
end

%% The solution
% So now what?  Clearly y_obs is not in the column space of A.  Yet no matter
% what we choose for b, we will get a set of predicted values that form a
% vector that is in the column-space.  Of course we want to find b so that
% the vector we get out is the closest one to the real solution.  Let's call
% the solution, ypred.  After all, there is no vector y, such that A*b =
% ypred.

% Now from the picture in front of you, it ought to be obvious that the
% closest we can get to y_obs is a vector that ends in the projection of
% y_obs
% on the blue plane.  I'm going to drop the end of vector y_obs onto the
% plane.  Try to ignore the math for the moment and grasp the geometry.

yproj = A * (inv(A' * A) * A' * y_obs)
hp = line([0 yproj(1)], [0 yproj(2)], [0 yproj(3)])
set(hp,'Color','g','LineStyle','--','LineWidth',3)


% The thick broken line is the projection of y_obs onto the column-space of
% A.  Let's call this vector yproj.  It is the predicted value that we would
% get out of A*b.  One way to see this is to return to the first plot and
% look at where bproj lies. These are the green asterisks connected by the
% green broken line.

figure(1)
hold off, plot(x,y_obs,'o'), hold on
plot(x,yproj,'g*', x,yproj,'g--')
set(gca,'XLim',[-3 3],'YLim',[-1 4])

% Our yproj is the column vector containing the ordinate values of these 
% three points.  You can see that these are the predicted values that lie
% on the best fitting line.  The way that we would say that this is the best 
% fitting line is that the distance from the asterisks to the data values 
% (vertically) is minimized.

% This distance is the difference between y_obs and A*b (i.e., ypred).
% Let's return to the vector diagram and convince ourselves.

%% Visualizing this error vector is the key to understanding least squares.

figure(2)

% It helps to construct a line that represents the error between b and
% bproj.  The error, E, is

E = y_obs - yproj

% Rather than displaying this vector emanating from the origin, we
% will add it to the end of bproj

hE = line([yproj(1) y_obs(1)], [yproj(2) y_obs(2)], [yproj(3) y_obs(3)])
set(hE,'Color','y','LineStyle','--','LineWidth',3)

% Do you appreciate that this E vector is orthogonal to the plane?
% There are lots of ways to talk about this, but you have to see
% and believe it first.  Don't go on unless you are absolutely comfortable
% with the idea that E is orthogonal to the plane spanned by the columns of 
% A.

% O.K., you made it.  Do you see that regardless of the dimension of the 
% problem, the least squares solution will be a projection of some
% high dimensional vector, y_obs, on some lower dimensional column space of A.  
% Sticking with the 2-dimensional case, if we were fitting a
% line to 4 points instead of 3, we would be projecting a 4-dimensional 
% vector, y_obs, onto the 2-dimensional plane spanned by the 2 4-dimensional
% column vectors of A.  One of these would be [1 1 1 1]'.  

%% Let's return to the matter of E and how to talk about it.  
% We have already said that it is orthogonal to the columns of A.
% That means that the dot product of E with any column of A should be 
% zero.  Since E is a column vector, we either have to take its transpose 
% and multiply E'*A(:,1) 

E' * A(:,1)     % 1st column 
E' * A(:,2)			% 2nd column

% or take the transpose of the column vectors and multiply these times
% the E column vector.

A(:,1)' * E		  % 1st column 
A(:,2)' * E			% 2nd column

% I'm assuming you know about dot products.  If you don't, you'll have to
% pick up a book. You might just recall that it is the product of the
% lengths times the cosine of the angle between the vectors.  So z'*z is the
% square of the length of z.  You can easily verify that the dot product of
% two orthogonal vectors is 0 e.g., visualize and compute

[1 0] * [0 1]'
[1 1] * [-1 1]'
[1 1] * [-2.7 2.7]'

% ... and so forth.
 
% If you don't understand dot products (sometimes called inner products)
% don't go on. 

% If you understand the notion of inner products and the fact that E is
% orthogonal to the columnspace of A, then we can move on.
% The first thing to say is that we can combine the dot products above
% into a single matrix operation.  Starting with E transpose, 

E' * A		

% should be a row vector with 2 0's. (You get very tiny numbers)
% Or we can transpose A and multiply it times the E column vector

A' * E		

% makes a column vector with 2 zeros. (within floating point precision)


% There is a very imporant idea here.  The Error vector, which is, after all,
% the residuals between the observed data and the fit, is orthogonal to
% the columnspace of A.  It lies in the left nullspace of A.  The reason 
% this is called the left nullspace is because we have to left-multiply 
% E * A to get our zeros.

% Once we see that this is so, we can derive the NORMAL EQUATIONS in matrix
% form.  We are looking for the best approximation to b in the equation A*b
% = y_obs.  We have already admitted that there is no vector b that will work
% because A*b always lies in the columnspace of A, whereas y_obs does not.
% So, instead, we'll solve for yproj So we say that what we really want is
% to minimize the error between our fit, A*b, and the observed values, y_obs.
% In least squares, we want to minimize the sum of the (y_obs - A*b).^2 What
% you should be absolutely convinced of is that y_obs - A*b is just the
% vector E, and that this vector is perpendicular to the columnspace of A.


% Start with E.  We have already convinced ourselves of the fact that A' * E
% = 0. If you don't remember this, go back to the last statement you
% executed.  But E is just 

% y_obs - A*b 

% Substituting for E, we get 

% A' * (y_obs - A*b) = 0 

% rearranging, we get 

% A'* y_obs = A'*A*b 

% Remember, we are trying to solve for b. So we multiply both sides of
% the equation by the inverse of (A'*A)

%		inv(A'*A) * A'* y_obs = inv(A'*A) * (A'*A) * b
%		inv(A'*A) * A'* y_obs = b

b = inv(A'*A) * A'* y_obs 	
% which is the normal equation (see Strang, p. 156, noting that our b is
% his xbar and our y_obs is his b)

% The column vector describes the best fitting line to the data. 
% 1.2857 is the intercept and 0.5714 is the slope. 
 
figure(1)
fplot('1.2857 + 0.5714*x',[-3 3],'c')

% Now there is one thing left to appreciate before we extend this to
% multiple regression. I promised that I would explain the projection matrix
% we used to make our pretty picture.  It would be nice to be able to
% project a vector onto the columnspace of A.  Well, we already have one
% example: A*b is the projection of y_obs onto the columnspace.  Now look at
% what you get when you left-multiply the equation for b by A.

% A * b = A * inv(A'*A) * A' * y_obs

% This leads to a simple idea.  The seemingly messy string of matrices,
% A * inv(A'*A) * A', projects y_obs onto the plane in our figure.

P = A * inv(A'*A) * A'

% Earlier in the tutorial, we used P to project the vertices of a square 
% onto the columnspace of A.  That's how we made the blue plane.

%% Extension to higher order polynomials
% So why is this so cool?  Because it gives us a very easy way to do
% multiple regression.  Suppose we have the ordered pairs,

data = [
     0     1
     1    53
     2    35
     3    -7
     4    87
     5   106
     6   250
     7   301
     8   346
     9   337
    10   256
    11   -25
    12  -421
]
		


figure(1)
hold off, plot(data(:,1),data(:,2),'o'),hold on

% We want to fit a 4th order polynomial.
% The A matrix has 5 columns

x = data(:,1);
y_obs = data(:,2);
A = [ones(size(x)) x x.^2 x.^3 x.^4]

% We want the best solution for xbar in the equation A*xbar = b

b = inv(A'*A) * A'* y_obs

% It's that simple.  We can get the fit back with
predvals = A*b
plot(x,predvals,'r*')

% and we can get a predicted curve by increasing the sampling density
xint = linspace(min(x),max(x),100)';
A = [ones(size(xint)) xint xint.^2 xint.^3 xint.^4];
predvals = A*b;
plot(xint,predvals,'r-')


%% If you want to go on, the next topic is weighted least squares
% Then, multiple regression with nondiagonal covariance matrix.
