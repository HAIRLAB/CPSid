function [wnewm2, p] = slrvar_update_weight_comp_m2(wold, ttr, xtr, alpha, lam)


D = length(wold);
p = zeros(D,1);

a = xtr'*(ttr-0.5);
wnewm2 = wold;
C = xtr'*2*diag(lam);

for dd = 1 : D
    aa = zeros(D,1);
    aa(dd) = alpha(dd);
    %B_dd = xtr'*2*diag(lam)*xtr(:,dd) + aa;
    B_dd = C*xtr(:,dd) + aa;
    e = a(dd) - B_dd' * wnewm2;
    wnewm2(dd) = wnewm2(dd) + e/B_dd(dd); 
    p(dd) = B_dd(dd);
end

