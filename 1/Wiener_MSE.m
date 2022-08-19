function J = Wiener_MSE(R, w, p, sd2)

R_tf = ishermitian(R);
d = eig(R);
if ~all(d >= 0) || ~R_tf
    error('Error. R is not positive semidefinite')
end

[R_dim1, ~] = size(R);
[w_dim1, ~] = size(w);
[p_dim1, ~] = size(p);
if R_dim1 ~= w_dim1 || w_dim1 ~= p_dim1
    error('Error. The dimensions of input arguments are not suitable')
end

if sd2<0 || ~isreal(sd2)
    error('Error. sd2 is negative or complex-valued')
end

J = sd2 - w'*p - p'*w + w'*R*w;

end
                
