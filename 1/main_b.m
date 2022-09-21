
% covariance matrix
R = [1.1 0.5 0.1;...
     0.5 1.1 0.5;...
     0.1 0.5 1.1];

% cross-correlation vector
p = [0.5;...
    -0.4;...
    -0.2];

% sd2 = E[|d(n)|^2]
sd2 = 1;

Re_w0 = linspace(-5, 5, 101);
J = zeros(1, 101);
for i = 1:101
    w = [Re_w0(i)+1i; -0.5+1i; -1];
    J(i) = real(Wiener_filter(R, w, p, sd2));
end

Jmin = min(J);
min_x = find(J == Jmin);

% plot
txt = ['optimal point is at ' 10 'Re\{w0\} = ' num2str(Re_w0(min_x)) 10 'and Jmin is ' num2str(Jmin)];
text(-1, 10, txt, 'color', 'red', 'Fontsize', 10)
hold on
plot(Re_w0(min_x), Jmin, '*', 'color', 'r', 'markersize',10);
plot(Re_w0, J, 'color', 'blue');
xlabel('Re\{w0\}');
ylabel('MSE');
legend('Optimal Point','MSE curve');


