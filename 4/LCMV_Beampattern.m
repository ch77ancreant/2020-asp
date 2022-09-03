

load('matX.mat');

% Covariance matrix estimation
R_hat = (matX*matX')/1000;
% Discrete grid theta -----------
theta = -90:0.01:90;

B_LCMV = zeros(1,length(theta));

% a(theta_s) ----------------------------
a_s = exp(1i*[0:10]*pi*sind(3.25)).';
% a(theta_i) ----------------------------
a_i = exp(1i*[0:10]*pi*sind(18.57)).';

C = [a_s , a_i];
g = [1; 10^(-5)];
w_LCMV = R_hat\C*((C'*(R_hat\C))\g);

% LCMV beamformer
for ii = 1:length(theta)
    a = exp(1i*[0:10]*pi*sind(theta(ii))).';
    B_LCMV(ii) = w_LCMV'*a;
end

% plot
plot(theta, 20*log10(abs(B_LCMV)), 'linewidth', 2);
hold on
line([3.25 3.25],[-80 20],'LineStyle',':','LineWidth',1.3,'Color',[0.8500 0.3250 0.0980])
line([18.57 18.57],[-80 20],'LineStyle',':','LineWidth',1.3,'Color',[0.4940 0.1840 0.5560])
legend('beampattern','signal direction','interference direction','Location','southwest')
title('LCMV Beampattern');
xlabel('$\theta$ (rad)','interpreter','latex');
ylabel('$|B_{\theta}(\theta)|$ (dB)','interpreter','latex');
axis([-90 90 -80 20]);
set(gca,'ytick',[-80 -60 -40 -20 0 20]);
set(gca,'xtick',[-90 -45 0 45 90]);
set(gca,'xticklabel',{'-\pi/2','-\pi/4','0','\pi/4','\pi/2'});
grid on



