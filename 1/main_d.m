
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

Re_w0 = linspace(-2, 2, 101);
Re_w2 = linspace(-2, 2, 101);
[xx, yy] = meshgrid(Re_w0, Re_w2);
J = zeros(101, 101);
for ii = 1:101
    for jj = 1:101
        w = [Re_w0(jj); -0.7683; Re_w2(ii)];
        J(ii, jj) = real(Wiener_MSE(R, w, p, sd2));
    end
end

Jmin =  min(min(J));
[min_y, min_x] = find(J == Jmin);

% plot
txt = [10 10 10 10 10 'optimal point is at ' 10 '(Re\{w0\}, Re\{w2\})=(' num2str(Re_w0(min_x)) ', ' num2str(Re_w2(min_y)) ')' 10 'and  Jmin is ' num2str(Jmin)];
text(0.2, Re_w2(min_y), txt, 'color', 'red', 'Fontsize', 10)
hold on
plot(Re_w0(min_x), Re_w2(min_y),'*', 'color', 'r', 'markersize',10);
contour(xx,yy,J,[0.35 0.6 1 2 3 4 5],'ShowText','on');
xlabel('Re\{w0\}');
ylabel('Re\{w2\}');


