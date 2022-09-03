clear;clc;close all

load('data.mat');

[N, L] = size(Y_tilde);
M = size(F_n1_n, 1);

% Initial condition
x_hat_n_given_ynm1 = zeros(M, 1);
K_n_nm1 = eye(M);

% Store the state vector
x_hat_n_given_yn = zeros(M, L);

% F(n,n+1) = inv(F(n+1,n))
F_n_n1 = inv(F_n1_n);

 for ii = 1:L
     
     % Kalman Gain
     R = C_n*K_n_nm1*C_n' + Q2_n;
     G = F_n1_n*K_n_nm1*C_n'*inv(R);
     
     % Riccati Equation
     K_n = K_n_nm1 - F_n_n1*G*C_n*K_n_nm1;
     K_n_nm1 = F_n1_n*K_n*F_n1_n' + Q1_n;
     
     % One-Step Predictor
     a = Y_tilde(:, ii) - C_n*x_hat_n_given_ynm1;
     x_hat_n_given_ynm1 = F_n1_n*x_hat_n_given_ynm1 + G*a;
     
     % Estimate the state vector at time n
     x_hat_n_given_yn(:, ii) = F_n_n1*x_hat_n_given_ynm1;
 end
 
 % plot
for k = 1:4
    figure(k);
    for m = 1:M
        subplot(M,1,m);
        if k == 1
            plot(1:L, real(x_hat_n_given_yn(m,:)));
            ylabel(['$Re\{\hat{x}_', num2str(m), '(n|\mathcal{Y}_n)\}$'],'interpreter','latex');
        elseif k == 2
            plot(1:L, imag(x_hat_n_given_yn(m,:)));
            ylabel(['$Im\{\hat{x}_', num2str(m), '(n|\mathcal{Y}_n)\}$'],'interpreter','latex');
        elseif k == 3
            plot(1:L, abs(x_hat_n_given_yn(m,:)));
            ylabel('magnitude');
        elseif k == 4
            plot(1:L, unwrap(angle(x_hat_n_given_yn(m,:))));
            ylabel('unwrapped phase','interpreter','latex');
        end
        title(['$\hat{x}_', num2str(m),'(n|\mathcal{Y}_n)$'],'interpreter','latex');
        xlabel('n');
        axis([0,L,-inf,inf]);
    end
end

%% Problem(c)
slope = zeros(M, 1);
for mm = 1:M
    A = [ones(L,1), (1:L).'];
    p = unwrap(angle(x_hat_n_given_yn(mm,:).'));
    a = (A.'*A)\A.'*p;
    slope(mm) = a(2);
end


 
 
 
 
 