C = [1 .2 .2; .2 1 .2; .2 .2 1]
Q = sqrtm(C)

n = 10000
x = normrnd(0,1,3,n);
y = Q*x;
size(y)
Cobs = cov(y')

[v lam] = eig(Cobs)




