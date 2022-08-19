
% covariance matrix
R = [1.1 0.5 0.1;...
     0.5 1.1 0.5;...
     0.1 0.5 1.1];

% cross-correlation vector
p = [0.5;...
    -0.4;...
    -0.2];

Re_w0 = linspace(-5, 5, 101);
J = zeros(1, 101);
for ii = 1:101
    w = [Re_w0(ii)+1i; -0.5+1i; -1];
    J(ii) = real(Wiener_MSE(R, w, p, sd2));
end

Jmin = min(J);
Jmin_x = find(J == Jmin);

% plot
txt = ['optimal point is at ' 10 'Re\{w0\} = ' num2str(Re_w0(Jmin_x)) 10 'and associated Jmin = ' num2str(Jmin)];
text(-1, 10, txt, 'color', 'red', 'Fontsize', 10)
hold on
plot(Re_w0(Jmin_x), Jmin, '*', 'color', 'r', 'markersize',10);
plot(Re_w0, J, 'color', 'blue');
xlabel('Re\{w0\}');
ylabel('MSE');
legend('Optimal Point','MSE curve')


