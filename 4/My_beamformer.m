clear; close all;

load('matX.mat');

% Covariance matrix estimation
R_hat = (matX*matX') / 1000;
a_s = exp(1i*[0:10]*pi*sind(3.25)).';
a_i = exp(1i*[0:10]*pi*sind(18.57)).';

% Update weights for every time
% and see output
w_d = zeros(11, 1000);
y_hat = zeros(1000, 1);

for k = 1:1000
    % Find s(t) and i(t) every time
    % C*y = x, y = inv(C'*C)*C'*x
    C = [a_s, a_i];
    x = matX(:, k);
    y = (C'*C)\C'*x;
    s_s = y(1);
    s_i = y(2);
    % Find noise
    n = x - a_s*s_s - a_i*s_i;
    
    % LCMV 
    D = [a_s, a_i, n];
    g = [1; 10^(-5); 10^(-5)];
    w_d(:,k) = R_hat\D*((D'*(R_hat\D))\g);
    y_hat(k) = w_d(:,k)'*x;  
end

% compare with MVDR beamformer
w_MVDR = (R_hat\a_s) / (a_s'*(R_hat\a_s));
y_MVDR = zeros(1, 1000);
for k = 1:1000
    y_MVDR(k) = w_MVDR' * matX(:,k);
end

% plot
subplot(2,1,1)
plot(1:1000, real(y_hat))
hold on
plot(1:1000, real(y_MVDR))
title('Beamformer Output (real part)')
ylabel('$Re\{\hat{y}(t)\}$', 'interpreter', 'latex')
xlabel('t','interpreter','latex')
legend('My beamformer','MVDR beamformer')
grid on

subplot(2,1,2)
plot(1:1000, imag(y_hat))
hold on
plot(1:1000, imag(y_MVDR))
title('Beamformer Output (imaginary part)')
ylabel('$Im\{\hat{y}(t)\}$', 'interpreter', 'latex')
xlabel('t','interpreter','latex')
legend('My beamformer','MVDR beamformer')
grid on


