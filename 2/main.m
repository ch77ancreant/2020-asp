
load('data_v.mat');

% H(z) -------------
d = [1, 1];
c = [1, -1/6, -1/6];

% x(n) = conv(v(n),h(n))
x = filter(d, c, v);

% compute the autocorrelation of x(n) by hand
% we can get rx(k) = 108/35*(1/2)^k - 18/35*(-1/3)^k
for k = 0:10
    r(k+1, 1) = 108/35*(1/2)^k - 18/35*(-1/3)^k;
end

[~, P, kappa] = Levinson_Durbin(r);

% x(n-1)
x_1 = [0, x(1:999)];

% stage 1, f1(n) and b1(n)
f1 = x + kappa(1)' * x_1;
b1 = x_1 + kappa(1) * x;

% Lattice Structure
% f1(n) ~ f10(n) in row1 ~ row10 of matrix f
% b1(n) ~ b10(n) in row1 ~ row10 of matrix b
f(1, :) = f1;
b(1, :) = b1;

for i = 2:10
    f(i, :) = f(i-1, :) + kappa(i)' * [0, b(i-1, 1:999)];
    b(i, :) = [0, b(i-1, 1:999)] + kappa(i) * f(i-1, :);
end

% Compute the prediction error power
Pf = zeros(1, 10);
Pb = zeros(1, 10);
for m = 1:10
    Pf(m) = sum(abs(f(m, :)).^2) / 1000;
    Pb(m) = sum(abs(b(m, :)).^2) / 1000;
end

% func is the PSD of x(n) --- get by hand
% integrate log(S(exp(j2pif))) and compute the prediction error bound
func = @(u) log10((2+2*cos(2*pi.*u)) ./ (19/18 - 5/18*cos(2*pi.*u) - 1/3*cos(4*pi.*u)));
Pmbound = exp(integral(func, -0.4999, 0.4999));

%plot
plot(1:10, Pf, '-*');
hold on 
plot(1:10, Pb, '-o');
plot(1:10, Pmbound*ones(1,10), 'r:', 'linewidth', 1.3);
axis([1 10 0.8 1.5]);
legend({'$\hat{P}_{f,m}$','$\hat{P}_{b,m}$','prediction error bound'},'interpreter','latex');
xlabel('m','interpreter','latex');
ylabel('prediction error power','interpreter','latex');




