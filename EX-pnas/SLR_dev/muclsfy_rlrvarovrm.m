function [ww, ix_eff_all, errTable_tr, errTable_te, parm, AXall, Ptr, Pte] =...
    muclsfy_rlrvarovrm(x_train, t_train, x_test, t_test, varargin)
% Multi-class classification by regularized logistic regression one-versus-rest classifiers.
% 
% Normalization, parameter estimation, and performance evaluation are executed.
%
% -- Usage
% [ww, ix_eff_all, errTable_Tr, errTable_te, parm, AXall, Ptr, Pte] =...
%   muclsfy_rlrvarovrm(x_train, t_train, x_test, t_test,  varargin)
%
% --- Input
% x_train :   [Nsamp_tr , Nfeat] 
% t_train :   [Nsamp_tr , 1]
% x_test  :   [Nsamp_te , Nfeat]
% t_test  :   [Nsamp_te , Nfeat]
%
% --- Optional Input
% parm = finputcheck(varargin, ...
%     {'scale_mode', 'string', {'all','each','stdall','stdeach','none'}, 'each';...
%      'mean_mode',  'string', {'all','each','none'}, 'each';...
%      'ax0',        'real',     [],  1;...
%      'nlearn',     'integer',  [1 inf],  1000;...
%      'nstep',      'integer',  [1 inf],  100;...
%      'amax',       'real',     [0 inf],  1e8;...
%      'usebias',    'boolean',  []     , 1;...
%      'norm_sep'  , 'boolean',  []     , 0;... 
%      'reduceparm', 'boolean',  []     , 0;...
%      'displaytext','boolean',  []     , 1;... 
%      'invhessian' ,'boolean',  []     , 0;...  
%      });
%
% --- Output
% ww          :   Estimated weight parameters. [Nfeat, Nclass]
% ix_eff_all  :   Index of features survived. cell array of {Nclass}
% errTable_tr :   Counting table of each label estimated. [Nclass, Nclass]
% errTbale_te :   Counting table of each label estimated. [Nclass, Nclass]
% parm        :   Parmaters used in this routine. [struct]
% AXall       :   History of hyperparameters updating. [Nfeat*Nclass Nlearn]
% Ptr         :   Probaility of observing every label in training data. [Nsamp_tr Nclass]
%                 This value is used to put a label on each sample.
% Pte         :   Probaility of observing every label in training data. [Nsamp_te Nclass]
%                 This value is used to put a label on each sample.
%
%
% 2016/03/30 OY
% * Bug fix of error table computation when t_test or t_train has single value
% 2009/08/10 OY
% * Bug fix when 't_test' has only a single label of 't_train'.
% 2009/06/01 OY
% * The function is renamed (previously run_rlr_bi_reg.m)
%
% 2007/12/28 OY
% * RLR version.
% 2006/10/23 OY
% * 'Nclass' is removed from inputs.
% 2006/09/12 OY
%  * A field "nprobe" is introduced.
% 2006/09/06 OY  
%  * Output format modified (error table as output)
%  * Comment modified
% 2006/08/02 OY
% 2006/05/26 OY
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.

if nargin < 4
    help muclsfy_rlrvarovrm;
    return
end

% char label -> number 
[t_train, label_names, Nclass] = label2num(t_train);
[t_test] = label2num(t_test, label_names);

[Nsamp_tr, Nfeat] = size(x_train);
Nsamp_te = size(x_test,1);

%% input check for optional parameter.
parm = finputcheck(varargin, ...
    {'scale_mode', 'string', {'all','each','stdall','stdeach','none'}, 'each';...
     'mean_mode',  'string', {'all','each','none'}, 'each';...
     'ax0',        'real',     [],  1;...
     'nlearn',     'integer',  [1 inf],  1000;...
     'nstep',      'integer',  [1 inf],  100;...
     'amax',       'real',     [0 inf],  1e8;...
     'usebias',    'boolean',  []     , 1;...
     'norm_sep'  , 'boolean',  []     , 0;... 
     'reduceparm', 'boolean',  []     , 0;...
     'displaytext','boolean',  []     , 1;... 
     'invhessian' ,'boolean',  []     , 0;...  
     });
 
if ~isstruct(parm)
   error(parm);
end

AMAX        = parm.amax;
ax0         = parm.ax0; 
Nlearn      = parm.nlearn;
Nstep       = parm.nstep;
usebias     = parm.usebias;
norm_sep    = parm.norm_sep;
displaytext = parm.displaytext;
invhessian  = parm.invhessian;
reduceparm  = parm.reduceparm;

%
if displaytext
    fprintf('-----------------------------------------------------\n');
    fprintf('  Multi-class classification by RLR-VAR one-vs-rest  \n');
    fprintf('-----------------------------------------------------\n');
end
% add a regressor for bias 
if usebias == 1
    Nfeat = Nfeat+1;
end
if reduceparm
    Nparm = (Nclass-1)*Nfeat;
else
    Nparm = Nclass*Nfeat;
end
% set ax0 
parm.ax0 = ax0;

% keep constant parameters
parm.nclass = Nclass;
parm.nparm = Nparm;
parm.nsamp_tr = Nsamp_tr;
parm.nsamp_te = Nsamp_te;
parm.nfeat    = Nfeat;

% normalize (sacling and baseline addjustment)
if norm_sep == 0
    [x_train, scale, base] = normalize_feature(x_train, parm.scale_mode, parm.mean_mode);
    [x_test, scale, base] = normalize_feature(x_test, parm.scale_mode, parm.mean_mode, scale, base);
else
    [x_train, scale, base] = normalize_feature(x_train, parm.scale_mode, parm.mean_mode);
    [x_test, scale, base] = normalize_feature(x_test, parm.scale_mode, parm.mean_mode);
end
% add a regressor for bias term
if usebias
    Xtr = [x_train, ones(Nsamp_tr,1)];
    Xte = [x_test, ones(Nsamp_te,1)];
else
    Xtr = [x_train];
    Xte = [x_test];
end

%---------------------
% VB
%---------------------
AXall = [];
ww = zeros(Nfeat,Nclass);
ix_eff_all = cell(Nclass,1);

for c = 1 : Nclass

    if displaytext, fprintf('\n Learning parameters for class %d .... \n', c); end

    if c == Nclass & reduceparm
        if displaytext, fprintf('All parameters fixed to zeros...\n'); end
        ww(:,c) = zeros(Nfeat,1);
        ix_eff_all{c} = [];
    else
        label = zeros(Nsamp_tr, 1);
        label(t_train == c) = 1;   % 1 for class c, 0 otherwise

        [w_e, ix_eff, W, AX] = rlr_learning_var2(label, Xtr,...
            'nlearn', Nlearn, 'nstep', Nstep, 'amax', AMAX, 'ax0', ax0, ...
            'invhessian', invhessian);

        ix_eff_all{c} = ix_eff;
        ww(:, c) = w_e;

        AXall = [AXall;AX];
    end
end

%-----------------------
% Training Correct
%-----------------------
[t_train_est, Ptr] = calc_label(Xtr, ww);

%-----------------------
% Test
%----------------------
[t_test_te, Pte] = calc_label(Xte, ww);

% remove baseline parameters from effective index
if usebias
    for cc = 1 : Nclass
        ix_eff_all{cc} = setdiff(ix_eff_all{cc}, Nfeat);
    end
end

%
errTable_tr = slr_error_table(t_train, t_train_est, unique([t_train;t_test]));
errTable_te = slr_error_table(t_test, t_test_est, unique([t_train;t_test]));

Pcorrect_tr = calc_percor(errTable_tr);
Pcorrect_te = calc_percor(errTable_te);

if displaytext,
fprintf(' Training Correct : %2.2f %%,  Test Correct : %2.2f %%\n', Pcorrect_tr, Pcorrect_te);
end