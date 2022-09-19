function J_n_RLS = ASP_RLS(x, d, lambda, delta, R)

J = zeros(R, size(d,2)-1);
J_n_RLS = zeros(1, size(d,2)-1);

for j = 1:R
    
    w_n = zeros(5, 1);
    x_n = [x(j, 1) 0 0 0 0].';
    P_n = delta^(-1).*eye(5);
    
    for i = 1:size(d,2)-1
        
        x_n(5) = []; xp = x_n;
        x_n = [x(j, i+1); xp];
        k_n = (lambda^(-1)*P_n*x_n)./(1 + lambda^(-1)*x_n'*P_n*x_n);
        Xi_n = d(j, i+1) - w_n'*x_n;
        J(j, i) = abs(Xi_n)^2;
        w_n = w_n + k_n*Xi_n';
        P_n = lambda^(-1)*P_n - lambda^(-1)*k_n*x_n'*P_n;
        
    end
end

for n = 1:size(d,2)-1
    J_n_RLS(n) = sum(J(:,n))/R; 
end

end

