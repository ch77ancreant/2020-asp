clear; close all;

load('matX.mat');

% Covariance matrix estimation
R_hat = (matX * matX') / 1000;

% Discrete grid theta -----------
theta = -90:0.01:90;

%% DOA estimation

% Scan all the theta and compute the spetrum
P_theta = zeros(1, length(theta));
for j = 1:length(theta)
    a = exp(1i*[0:10]*pi*sind(theta(j))).';
    P_theta(j) = 1/(a'*(R_hat\a));
end

% Find the peaks on the sepctrum
[pks, locs] = findpeaks(abs(P_theta));
locs = theta(locs);

% plot
figure
plot(theta, abs(P_theta))
title('MVDR spectrum')
xlabel('$\theta$ (rad)','interpreter','latex')
set(gca,'xtick',[-90 -45 0 45 90])
set(gca,'xticklabel',{'-\pi/2','-\pi/4','0','\pi/4','\pi/2'})
ylabel('$P_{MVDR}(\theta)$','interpreter','latex')
axis([-90 90 0 8])
grid on


%% Beampattern

B_Uniform = zeros(1, length(theta));
B_AS = zeros(1, length(theta));
B_MVDR = zeros(1, length(theta));

% a(theta_s) ----------------------------
a_s = exp(1i*[0:10]*pi*sind(3.25)).';

% -----------------------------------------
% The weights of different beamformer
w_Uniform = 1/11 * ones(11,1);
w_AS = 1/11 * a_s;
w_MVDR = (R_hat\a_s) / (a_s'*(R_hat\a_s));

% B_theta(theta) = w' x a(theta)
% Uniform weights, w = 1/N
for i = 1:length(theta)
     B_Uniform(i) =  sum(exp(1i*pi*sind(theta(i))*[0:10]));
end
 
% Array steering with source DOA
for i = 1:length(theta)
     B_AS(i) =  sum(exp(1i*pi*(sind(theta(i)) - sind(3.25))*[0:10]));
end
  
% MVDR beamformer with source DOA
for i = 1:length(theta)
    a = exp(1i*[0:10]*pi*sind(theta(i))).';
    B_MVDR(i) = w_MVDR'*a;
end

Beam_ULA = 1/11 * abs(B_Uniform);
Beam_AS = 1/11 * abs(B_AS);
Beam_MVDR = abs(B_MVDR);

% plot 
figure
plot(theta, 20*log10(Beam_ULA),'linewidth',2)
hold on
plot(theta, 20*log10(Beam_AS),'linewidth',2)
plot(theta, 20*log10(Beam_MVDR),'linewidth',2)
title('Beampattern')
xlabel('$\theta$ (rad)','interpreter','latex')
ylabel('$|B_{\theta}(\theta)|$ (dB)','interpreter','latex')
legend('Uniform weights','Array steering','MVDR')
axis([-90 90 -80 20])
set(gca,'ytick',[-80 -60 -40 -20 0 20])
set(gca,'xtick',[-90 -45 0 45 90])
set(gca,'xticklabel',{'-\pi/2','-\pi/4','0','\pi/4','\pi/2'})
grid on


%% Beamformer output

y_Uniform = zeros(1, 1000);
y_AS = zeros(1, 1000);
y_MVDR = zeros(1, 1000);

% Uniform weight
for k = 1:1000
    y_Uniform(k) = w_Uniform' * matX(:,k);
end

% Array steering
for k = 1:1000
    y_AS(k) = w_AS' * matX(:,k);
end

% MVDR beamformer
for k = 1:1000
    y_MVDR(k) = w_MVDR' * matX(:,k);
end

% plot
figure
subplot(2,1,1)
plot(1:1000, real(y_Uniform))
hold on
plot(1:1000, real(y_AS))
plot(1:1000, real(y_MVDR))
title('The output of beamformer (real parts)')
xlabel('t','interpreter','latex')
ylabel('$Re\{y(t)\}$','interpreter','latex')
legend('Uniform weights','Array steering','MVDR')

subplot(2,1,2)
plot(1:1000, imag(y_Uniform))
hold on
plot(1:1000, imag(y_AS))
plot(1:1000, imag(y_MVDR))
title('The output of beamformer (imaginary parts)')
xlabel('t','interpreter','latex')
ylabel('$Im\{y(t)\}$','interpreter','latex')
legend('Uniform weights','Array steering','MVDR')


