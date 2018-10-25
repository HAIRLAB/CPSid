function kappa = calc_kappa(tab),
%
% kappa 0 ---- 0.4 ---- 0.6 --- 0.8 --- 1
%         low      modest   high   perfect
% 2008 OY

N = sum(tab(:));
Ncor = sum(diag(tab));

Njud1 = sum(tab,2)';
Njud2 = sum(tab,1);

pcor_coincidence = sum(Njud1/N .* Njud2/N) ;
pcor_judge = Ncor/N;

kappa = (pcor_judge - pcor_coincidence)/(1-pcor_coincidence);
