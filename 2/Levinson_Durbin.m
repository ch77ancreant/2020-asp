function [a, P, kappa] = Levinson_Durbin(r)

% initial value
P0 = r(1); 
d0 = conj(r(2));
%-----------------------------
kappa = zeros(size(r,1)-1, 1);
Pp = zeros(size(r,1)-1, 1);
d = zeros(size(r,1)-2, 1);
%-----------------------------
kappa(1) = -d0/P0;
Pp(1) = P0*(1-abs(kappa(1))^2);
d(1) = conj(r(3)) + kappa(1) * conj(r(2));
a{1} = [1; kappa(1)];
%-----------------------------
% Levinson-Durbin Algorithm
for ii = 2:size(r,1)-1
    kappa(ii) = -d(ii-1)/Pp(ii-1);
    Pp(ii) = Pp(ii-1)*(1-abs(kappa(ii))^2);
    ak = a{ii-1};
    a{ii} = vertcat(ak,0) + kappa(ii) * vertcat(0,flipud(conj(ak)));
    al = a{ii};
    if ii<=size(r,1)-2
        d(ii) = fliplr(r(2:(ii+2))')*al;
    end
end
P = [P0;Pp];


