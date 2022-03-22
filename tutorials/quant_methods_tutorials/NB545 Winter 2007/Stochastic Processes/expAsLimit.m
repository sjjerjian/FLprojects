clear all

k = .1
m = [1:100]';
y = (1 - k./m).^m;

figure(10), clf
plot(m,y,'k-','LineWidth',3)
hold on
line(get(gca,'XLim'), exp(-k)*[1 1],'Color','r')
axprefs(gca)
xlabel('m')
ylabel('y')


%% look at the derivative wrt k
dx = .01
xmax = 10;
x = [0:dx:xmax]';
m = 1000
k = -1
y = (1 - k*x./m).^m;
dydx = diff(y)./dx;
figure(10), clf
plot(x,y,'k',x,exp(-k*x),'y--',x(2:end),dydx,'r',x,-k*exp(-k*x),'g--')

%% Use symbolic math toolbox
syms x n real
f = (1+(x/n))^n
fdot = diff(f,'x')
pretty(fdot)
% anyway it's obvious
% deriv of f is (1 + (x/n))^(n-1)
% = (1 + (x/n))^n / (1 + (x/n))
% = f / (1 + (x/n))
% as n gets large, the denominatory --> 1
% fdot -> f
% which implies f -> exp(x)




