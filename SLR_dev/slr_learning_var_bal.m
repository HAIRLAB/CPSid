function [w, ix_eff, W, AX] = slr_learning_var_bal(label, X, varargin)
% Learning parameters of ARD-sparse logistic regression model with
% balanced weight to compensate the effect of unbalanced samples (SLR-var-bal).
%
% The estimation algorthm is derived from lower bound method. 
% The likelihood function is approximated by Gaussian distribution using variational parameters.  (SLR-var)  
%
% Note that label vectors must be consisting of {0,1} values. 
%
% -- Usage 
% [w, ix_eff, W, AX] = slr_learning_var_bal(label, X, varargin)
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
% nohesssian : If 0 (FALSE), the diagonal elements of posterior variance is calculated without   
%         doing matrix inversion of size #feature*#feature  
% balanced : If 1, each sample is weighted. In default setting, balanced_weight is specfied so that the weight is proportional to 
%          inverse of number of samples in each label.
% balanced_weight : Weigth for class 1 samples. Valid only when 'balanced = 1'. 
%
% -- Output
% w_e    : Estimated parameters
% ix_eff : Index of the effective feature
% W      : History of parameter learning
% AX     : History of hyper parameter learning
%
% -- Example
% > [w_e, ix_eff, W, AX] = slr_learning_var_bal(t, X,...
%   'nlearn', 100, 'nstep', 10);
%
% 2015/06/29 By O.Yamashita
% * modify weighting scheme. Adding normalization steps.
% 2012/04/20 By O.Yamashita
% * a varibale name is modified.  'nohessian' -> 'invhessian'
% 2008/12/18 by Okito Yamashita
% * Weigthing samples in case of unbalanced samples. 'balanced'
% 2007/01/10 By Okito Yamashita
% * A field name 'fast' is changed to 'nohesssian', and also its variable type.
% 2006/10/04 by Okito Yamashita 

%% Error check
if nargin < 2
    help slr_learning_var_bal
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
     'invhessian',  'boolean', [],      0;...
     'balanced',   'boolean', [],      0;...
     'balanced_weight', 'real', [0 inf], NaN;...
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
invhessian   = opt.invhessian;
balanced = opt.balanced;
balanced_weight = opt.balanced_weight;

% In case of unbalanced samples, weighting samples is considered.
if ~balanced,
    v = zeros(Nsamp,1);
    ix0 = find(label == 0);
    ix1 = find(label == 1);
    if isnan(balanced_weight)   
        % Default weights : inverse of # of samples
        v(ix0) = 1/(length(ix0));
        v(ix1) = 1/(length(ix1));
    else
        v(ix0) = 1;
        v(ix1) = balanced_weight; 
    end
    v = v / (sum(v)) * Nsamp;  % normalize weights so that sum of weights become the total number of samples. 
else
    v = ones(Nsamp,1);
end

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

    
    if invhessian  
        %-----------------------------------------------
        % Calculate Hessian directly. 
        % Require inverse of matrix size (Nfeat*Nfeat)
        % This mode should be fast when Nfeat << Nsamp
        %-----------------------------------------------
        
        %H_eff = diag(ax_eff) + 2*X_eff'*diag(lam)*X_eff;
        lamv = lam .* v;
        H_eff = diag(ax_eff) + 2*X_eff'* (lamv(:,ones(1,length(ix_eff))).*X_eff);
        S_eff = inv(H_eff);
        % W-step
        w_eff = 1/2*S_eff*X_eff'*(v.*(2*label-1));
        
        %% A-step
        ax_eff = (1-ax_eff.*diag(S_eff))./(w_eff.^2);
    
      %  ax_eff = 1 ./ ((w_eff.^2)+diag(S_eff));
        %% Xi-step
        Xw = X_eff*w_eff;
        xi2 = Xw.^2 + sum(X_eff.*(X_eff*S_eff),2);
    else %%%%% New implementation  
        %-----------------------------------------------
        % Does not compute Hessian directly. 
        % Require inverse of matrix size (Nsamp*Nsamp)
        % This mode should be fast when Nfeat >> Nsamp
        %-----------------------------------------------
        ia = 1./ax_eff;
        ilam = 1./(lam.*v);
        iAX = ia(:,ones(1,Nsamp)) .* X_eff';
        XiAX = X_eff*iAX;

        C = 1/2*diag(ilam) + XiAX;
        iC = inv(C);

        ib = zeros(length(ix_eff),1);
        for jj = 1 : length(ix_eff)
            ib(jj) = iAX(jj,:)*iC*iAX(jj,:)';
        end

        dS_eff = ia - ib;
        ilamy =   ilam.*(v.*(2*label-1));

        % W-step
        w_eff = 1/4* iAX * iC* ilamy;

        %% A-step
        ax_eff = (1-ax_eff.*dS_eff)./(w_eff.^2);
        %% Xi-step
        Xw = X_eff*w_eff;

        xi2 = Xw.^2 + 1/2*diag(XiAX*iC*diag(ilam));
        
    end
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
       display('Caution: No feature is survived !!');
       break;
   end
   
end

