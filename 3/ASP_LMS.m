function J_n_LMS = ASP_LMS(x, d, Mu, R)

J = zeros(R, size(d,2));
J_n_LMS = zeros(1, size(d,2));

for j = 1:R
    
    w_n = zeros(5, 1);
    x_n = [x(j, 1) 0 0 0 0].';
    e_n = d(j, 1);
    J(j, 1) = abs(e_n)^2;

    for i = 2:size(d,2)
        
        w_n = w_n + Mu*x_n*e_n';
        x_n(5) = []; xp = x_n;
        x_n = [x(j, i); xp];
        e_n = d(j, i) - w_n'*x_n;
        J(j, i) = abs(e_n)^2;
        
    end
end

for n = 1:size(d,2)
    J_n_LMS(n) = sum(J(:,n))/R; 
end

end

