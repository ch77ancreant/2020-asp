
% x(n): complex input process
% R: covariance matrix. R = E[x(n)x^H(n)]
R = [1.1 0.5 0.1;...
     0.5 1.1 0.5;...
     0.1 0.5 1.1];

% p: cross-correlation vector. p = E[x(n)d*(n)]
p = [0.5;...
    -0.4;...
    -0.2];

% d(n): complex desired process. 
% sd2 = E[|d(n)|^2]
sd2 = 1;

% weight vector (optimal)
w_opt = R \ p;

% ---------------------------------
% output, MSE J
J = Wiener_MSE(R, w_opt, p, sd2)

