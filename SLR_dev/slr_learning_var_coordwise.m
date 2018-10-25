function [w, ix_eff, W, AX] = slr_learning_var_coordwise(label, X, varargin)
% Learning parameters of ARD-sparse logistic regression model with the
% component-wise approximation (SLR-var-comp).
%
% The estimation algorthm is derived from lower bound method. 
% The likelihood function is approximated by Gaussian distribution using variational parameters.  (SLR-var)  
%
% Note that label vectors must be consisting of {0,1} values. 
%
% -- Usage 
% [w, ix_eff, W, AX] = slr_learning_var_coordwise(label, X, varargin)
% 
% -- Input
% label : Teacher label vector consisting of {0,1}  [#sample*1]
% X     : Explanatory matrix                        [#sample*#feature]
%
% -- Field of Optional Input
% ax0   : Initial value of relevance parameters ax
% xi0   : Initial value of variational parameters xi
% nlearn : # of learning 
% nstep  : # of step at which parameters updataing is kept.
% amax   : Truncation criteria. Parameters whose relevance paramater ax is l
%          arger than this value are eliminated from the further iterations.
%
% -- Output
% w_e    : Estimated parameters
% ix_eff : Index of the effective feature
% W      : History of parameter learning
% AX     : History of hyper parameter learning
%
%
% 
% 2009/08/19 OY
% * modify a comment when no feature is survived
% 2009/07/09 OY
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.

%% Error check
if nargin < 2
    help slr_learning_var_coordwise
end

%% # of parameters     
[Nsamp,Nparm]  = size(X); 

%% input check for optional parameter.
opt = finputcheck(varargin, ...
    {'ax0',         'real',    [],      ones(Nparm,1);...
     'xi0',         'real',    [],      2*ones(Nsamp,1);...
     'nlearn',      'integer', [1 inf], 150;...
     'nstep',       'integer', [1 inf], 10;...
     'amax',        'real',    [0 inf], 1e8;...
     });
 
if ~isstruct(opt)
   error(opt);
end
    
if length(opt.ax0) ~= Nparm 
    error(['ax0 must be a vector of size ' num2str(Nparm), 'x 1 !!']);
end

if length(opt.xi0) ~= Nsamp 
    error(['xi0 must be a vector of size ' num2str(Nsamp), 'x 1 !!']);
end


%% Initial value for A-step and W-step
Nlearn = opt.nlearn;
Nstep  = opt.nstep;
AMAX   = opt.amax;

W = [];
AX = [];

ax = opt.ax0;   
xi = opt.xi0;
w  = zeros(Nparm,1); 
ix_eff = [1:Nparm]'; %effective index 


for nlearning = 1 : Nlearn
    
    %% Effective parameters
    ax_eff = ax(ix_eff);
    w_eff = w(ix_eff);
    X_eff = X(:, ix_eff);
    
    %% W-step
    lam = tanh(xi/2)./xi/4;  

    % w_eff, p_eff : mean and precision of each weight
    [w_eff, p_eff] = slrvar_update_weight_comp_m2(w_eff, label, X_eff, ax_eff, lam);
    
    
    %% A-step
    ax_eff = (1-ax_eff./p_eff)./(w_eff.^2);

    %  ax_eff = 1 ./ ((w_eff.^2)+1./p_eff);
        
    %% Xi-step
    Xw = X_eff*w_eff;
    xi2 = Xw.^2 + sum(X_eff.*(X_eff*diag(1./p_eff)),2);
    
    
    
    %   %%%%%%%%%%%%%%%%%%%%
    %
    %   subplot(3,1,1),
    %   plot(dS_eff - diag(S_eff));
    %   subplot(3,1,2),
    %   plot(w_eff - w_eff2);
    %
    %   subplot(3,1,3)
    %   plot(xi2 - xi22);
    %     pause(2)
    
  
    
    xi = sqrt(xi2);
     
    %% Prune ineffective parameters
    w = zeros(Nparm,1);
    w(ix_eff) = w_eff;
    ax(ix_eff) = ax_eff;
    ix_eff = find(ax < AMAX);
  
  
   %% Keep history of parameter updating
   
   if mod(nlearning, Nstep) == 0
        fprintf('Iterations : %d, Feature Remained: %d \n', nlearning, length(ix_eff)); 
        W(:, nlearning/Nstep) = w;
        AX(:,nlearning/Nstep) = ax;
        
     %   semilogy(ax);
     %   pause(0.5);
   end
   
   if isempty(ix_eff)
       display('Caution: No feature is survived !! The weight vector is set to zeros. ');
       w = zeros(Nparm,1);
       break;
   end
   
end

AX = [opt.ax0 AX];
