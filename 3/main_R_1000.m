clear;clc;close all;

load('matV.mat');
d = matV;

% Monte-Carlo runs
R = 1000;

% H(z) -----------------------------------
b = [1, 1/10];
a = [1, -1/6, -1/6];

% input signal 'x' and desired signal 'd'
x = zeros(1000, 500);
for k = 1:1000
    x(k,:) = filter(b, a, matV(k,:));
end


% implement ------------------------------
J_n_LMS1 = ASP_LMS(x, d, 0.1, R);
J_n_LMS2 = ASP_LMS(x, d, 0.2, R);
J_n_NLMS1 = ASP_NLMS(x, d, 0.2, R);
J_n_NLMS2 = ASP_NLMS(x, d, 0.8, R);
J_n_RLS1 = ASP_RLS(x, d, 0.75, 0.01, R);
J_n_RLS2 = ASP_RLS(x, d, 0.75, 0.1, R);
J_n_RLS3 = ASP_RLS(x, d, 0.95, 0.01, R);

%% plot -------------------------------------
semilogy(1:500, J_n_LMS1, 'Marker','.');
hold on
semilogy(1:500, J_n_LMS2, 'Marker','.');
semilogy(1:500, J_n_NLMS1, 'Marker','.');
semilogy(1:500, J_n_NLMS2, 'Marker','.');
semilogy(1:499, J_n_RLS1, 'Marker','.');
semilogy(1:499, J_n_RLS2, 'Marker','.');
semilogy(1:499, J_n_RLS3, 'Marker','.');
legend('LMS: \mu = 0.1','LMS: \mu = 0.2','NLMS: \mu = 0.2','NLMS: \mu = 0.8',...
    'RLS: \lambda = 0.75, \delta = 0.01','RLS: \lambda = 0.75, \delta = 0.1','RLS: \lambda = 0.95, \delta = 0.01');
xlabel('$n$','interpreter','latex');
ylabel('$\hat{J}(n)$','interpreter','latex');
title('Learning curves with R = 1000');
grid on
