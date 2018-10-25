% 2 class classfication by linear SLR
% Simple classification when multiple features are relevant but either of
% them is not so strongly relevant.
% 
% discriminant fucntion : linear function 
% training data : simulated data generated from Gaussian Mixtures 
% test data     : simulated data generated from Gaussian Mixtures 

clear
close all

data = 1;

fprintf('This is a demo how the sparse logistic regression model works...\n');
%----------------------------
% Generate Data
%----------------------------
switch data
     case 1,

        D = 400;
        Ntr = 200;
        Nte = 100;
        % mean
        mu1 = zeros(D,1);
        mu2 = [1.5; 0; zeros(D-2,1)];
        % covariance
        S = diag(ones(D,1));
        ro = 0.8;
        S(1,2) = ro;
        S(2,1) = ro;


        [ttr, xtr, tte, xte, g] = gen_simudata([mu1 mu2], S, Ntr, Nte);

        fprintf('\nThe data is generated from 2 Gaussian Mixture model of which centers (mean) are different.\n');
        fprintf('Input feature dimension is %d. \n',D);
        fprintf('But only the first dimension has difference in mean value between two classes,\n');
        fprintf('and the other dimension has same mean value.\n');
        fprintf('Therefore only the first dimension is detected as a meaningful feature,\n')
        fprintf('if you select features by feature-wise t-value ranking method.\n');
        fprintf('However due to the correlataion between the second dimension and the first dimension,\n')
        fprintf('inclusion of the second dimension makes classfication more accurate.\n');
        fprintf('For comparison, this demo also computes the classification performance of linear RVM, \n');
        fprintf('which is Bayesian counterpart of support vector machine (SVM).');

    case 2,
        D = 100;
        Ntr = 200;
        Nte = 100;
        % mean
        mu1 = zeros(D,1);
        mu2 = [[1:-0.02:0]'; zeros(D-51,1)];
        % covariance
        S = diag(ones(D,1));

        [ttr, xtr, tte, xte, g] = gen_simudata([mu1 mu2], S, Ntr, Nte);

        fprintf('\nThe data is generated from 2 Gaussian Mixture model of which centers (mean) are different.\n');
        fprintf('Input feature dimension is %d. \n',D);
        fprintf('The first 50 dimension has difference in mean value between two classes,\n');
        fprintf('and the other dimension has same mean value.\n');
        fprintf('The degree of relevance in the first 50 dimensions are manipulated by the mean values in class 1.\n')
        fprintf('For comparison, this demo also computes the classification performance of linear RVM, \n');
        fprintf('which is Bayesian counterpart of support vector machine (SVM).');
    
    case 3,
        load('../TESTDATA/real_binary', 'TRAIN_DATA', 'TEST_DATA', 'TRAIN_LABEL', 'TEST_LABEL');
        ttr = TRAIN_LABEL;
        tte = TEST_LABEL;
        xtr = TRAIN_DATA;
        xte = TEST_DATA;
        [Ntr,D] = size(TRAIN_DATA);
        [Nte] = size(TEST_DATA,1);
        
        fprintf('\nThe data is generated from a real experimental EEG data.\n');
        fprintf('In the experiment a subject executed either left or right finger tapping.\n');
        fprintf('The data has already been processed appropriately for classification.\n');
end

%--------------------------------
% Plot data (First 2 dimension)
%--------------------------------
slr_view_data(ttr, xtr);
axis equal;
title('Training Data')

% fprintf('\n\nPress any key to proceed \n\n');
% pause

%--------------------------------
% Learn Paramters
%--------------------------------

gamma = 15;

tic
fprintf('\nL1-Laplace-comp!!\n')
[ww_oc, ix_eff_oc, errTable_tr_oc, errTable_te_oc, parm,  ptr_oc, pte_oc] =...
    biclsfy_l1slrc(xtr, ttr, xte, tte, gamma,...
     'nlearn', 300, 'mean_mode', 'none', 'scale_mode', 'none');
toc

tic
fprintf('\nL1-Laplace!!\n')
[ww_o, ix_eff_o, errTable_tr_o, errTable_te_o, parm, ax_o, ptr_o, pte_o] =...
    biclsfy_l1slrlap(xtr, ttr, xte, tte, gamma^2,...
    'wdisp_mode', 'off', 'nlearn', 300, 'mean_mode', 'none', 'scale_mode', 'none');
toc

tic
fprintf('\nOLD version (ARD-Laplace)!!\n')
[ww_on, ix_eff_on, errTable_tr_on, errTable_te_on, parm_n, ax_on, ptr_on, pte_on] = biclsfy_slrlap(xtr, ttr, xte, tte,...
    'wdisp_mode', 'off', 'nlearn', 300, 'mean_mode', 'none', 'scale_mode', 'none');
toc

%% differenece

% ww_o(:,1) + ww_on
% ix_eff_o{1} - ix_eff_on{1}
% errTable_tr_o - errTable_tr_on
% errTable_te_o - errTable_te_on
% ax_o - ax_on
% ptr_o - ptr_on
% pte_o - pte_on

% tic
% fprintf('\nFast version (ARD-Variational)!!\n')
% [ww_f, ix_eff_f, errTable_tr_f, errTable_te_f, parm, AXall_f] = run_smlr_bi_var(xtr, ttr, xte, tte,...
%     'nlearn', 300, 'mean_mode', 'none', 'scale_mode', 'none');
% toc
% 
% tic
% fprintf('\nNew Fast version (ARD-Variational)!!\n')
% [ww_fn, ix_eff_fn, errTable_tr_fn, errTable_te_fn, parm, AXall_fn] = biclsfy_slrvar(xtr, ttr, xte, tte,...
%     'nlearn', 300, 'mean_mode', 'none', 'scale_mode', 'none', 'invhessian', 1);
% toc


% tic
% fprintf('\n\nLinear RVM (ARD-Variational)!!\n')
%     [ww_rvm, ix_eff_rvm, errTable_tr_rvm, errTable_te_rvm, g_rvm] = run_rvm(xtr, ttr, xte, tte, 0, ...
%         'nlearn', 1000, 'nstep', 1000, 'mean_mode', 'none', 'scale_mode', 'none', 'amax', 1e8);
% toc

