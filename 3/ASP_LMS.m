function J_n_LMS = ASP_LMS(x, d, Mu, R)

J = zeros(R, 500);
J_n_LMS = zeros(1, 500);

for jj = 1:R
    
    w_n = zeros(5, 1);
    x_n = [x(jj, 1) 0 0 0 0].';
    e_n = d(jj, 1);
    J(jj, 1) = abs(e_n)^2;

    for ii = 2:500
        
        w_n = w_n + Mu*x_n*e_n';
        x_n(5) = []; xp = x_n;
        x_n = [x(jj, ii); xp];
        e_n = d(jj, ii) - w_n'*x_n;
        J(jj, ii) = abs(e_n)^2;
        
    end
end

for n = 1:500
    J_n_LMS(n) = sum(J(:,n))/R; 
end

end

