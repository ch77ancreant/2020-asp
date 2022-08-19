
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
Im_w0 = linspace(-4 ,4, 81);
[xx, yy] = meshgrid(Re_w0, Im_w0);
J = zeros(81, 101);
for ii = 1:81
    for jj = 1:101
        w = [Re_w0(jj) + 1i*Im_w0(ii); 1 + 1i; 0.5];
        J(ii, jj) = real(Wiener_MSE(R, w, p, sd2));
    end
end

Jmin = min(min(J));
[min_y, min_x] = find(J == Jmin);
fprintf('optimal point is at (Re{w0}, Im{w0}) = (%.2f, %.2f)\n', Re_w0(min_x), Im_w0(min_y))
fprintf('and Jmin is %.2f\n', Jmin)

% plot
surf(xx, yy, J);
colorbar;
hold on
plot3(Im_w0(min_y), Re_w0(min_x),  Jmin, '*', 'color', 'r',...
        'Markersize', 10, 'Linewidth', 2);
xlabel('Re\{w0\}');
ylabel('Im\{w0\}');
zlabel('MSE');


     