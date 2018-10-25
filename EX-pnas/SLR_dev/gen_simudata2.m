function [ttr, xtr, tte, xte, g] = gen_simudata2(MU, S, Ntrs, Ntes)
% Generate simulation data2
%
% MU : 
% S  : 
% Ntrs : number of samples for each class [Nclass*1]
% Ntes : number of smaples for each class [Nclass*1]
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.
  
[Nfeat, Nclass] = size(MU);

if ndims(S) == 2,
    SS = repmat(S,[1,1,Nclass]);
else 
    SS = S;
end


if length(Ntrs) == 1,  NTRS = repmat(Ntrs, [Nclass,1]); else, NTRS = Ntrs; end
if length(Ntes) == 1,  NTES = repmat(Ntes, [Nclass,1]); else, NTES = Ntes; end


% simulation data generation
xtr = [];
ttr = [];
for ii = 1 : Nclass
   mutmp = MU(:,ii);
   Stmp  = squeeze(SS(:,:,ii));
   xtmp  = randmn(mutmp, Stmp, NTRS(ii));
   xtr = [xtr; xtmp'];
   ttr = [ttr; ii*ones(NTRS(ii),1)];
end

% simulation data generation
xte = [];
tte = [];
for ii = 1 : Nclass
   mutmp = MU(:,ii);
   Stmp  = squeeze(SS(:,:,ii));
   xtmp  = randmn(mutmp, Stmp, NTES(ii));
   xte = [xte; xtmp'];
   tte = [tte; ii*ones(NTES(ii),1)];
end
   


if Nclass == 2,
    ttr = ttr - 1;
    tte = tte - 1;
end

g.MU = MU;
g.S  = S;
g.Ntr = NTRS;
g.Nte = NTES;
g.Nclass = Nclass;
g.Nfeat  = Nfeat;
