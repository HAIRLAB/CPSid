% Multi class classfication demos
% 
% discriminant fucntion : linear function 

clear
close all

data = 1;  % 1 : artificial data 2 : real data

fprintf('This demo shows how the multi-class classification problems is solved ...\n');
%----------------------------
% Generate Data
%----------------------------
switch data
     case 1,
         N = 100; % sample per class
         
         fprintf('3 classes and 100 samples are used \n');

         mu = [-2 2 ; 0 0 ; 2 2 ;];
         sig = [1; 1; 1;];
       
         Nclass = length(sig);
         mark = {'ro', 'b+', 'g.'};

         X = [];
         label = [];
         for c = 1 : Nclass
             x = randmn(mu(c,:)', sig(c), N)';
             X = [X; x];
             label = [label; c*ones(N,1)];
         end

         [ixtr,ixte] = separate_train_test(label, 0.5);
         xtr = X(ixtr,:);
         xte = X(ixte,:);
         ttr = label(ixtr,:);
         tte = label(ixte,:);
         

    case 2,
        load('../TESTDATA/real_fourclass', 'TRAIN_DATA', 'TEST_DATA', 'TRAIN_LABEL', 'TEST_LABEL');
        ttr = TRAIN_LABEL;
        tte = TEST_LABEL;
        xtr = TRAIN_DATA;
        xte = TEST_DATA;
        [Ntr,D] = size(TRAIN_DATA);
        [Nte] = size(TEST_DATA,1);
        
        fprintf('\nThe data is generated from a real experimental EEG data.\n');
        fprintf('In the experiment a subject imagened either of left-hand, right hand, foot or tongue movement.\n')
        fprintf('The data has already been processed appropriately for classification.\n');

    otherwise,
        error('Choose data = 1 or data = 2 ... !');

end

%--------------------------------
% Plot data (First 2 dimension)
%--------------------------------
slr_view_data_multi(ttr, xtr);
axis equal;
title('Training Data')

fprintf('\n\nPress any key to proceed \n\n');
pause(0.1)

%--------------------------------
% Learn Paramters
%--------------------------------

tic
fprintf('\n\nSMLR (ARD-Laplace)!!\n')
    [ww_s, ix_eff_s, errTable_tr_s, errTable_te_s, g_s, ax, ptr, pte] = run_smlr(xtr, ttr, xte, tte, ...
        'nlearn', 300, 'nstep', 100, 'wdisp_mode', 'off', ...
        'mean_mode', 'none', 'scale_mode', 'none', 'amax', 1e8);
toc


tic
fprintf('\n\nSMLR new (ARD-Laplace)!!\n')
    [ww_sn, ix_eff_sn, errTable_tr_sn, errTable_te_sn, g_sn, axn, ptrn, pten] =...
        muclsfy_smlr(xtr, ttr, xte, tte, ...
        'nlearn', 300, 'nstep', 100, 'wdisp_mode', 'off', ...
        'mean_mode', 'none', 'scale_mode', 'none', 'amax', 1e8);
toc

ww_s - ww_sn
errTable_tr_s - errTable_tr_sn
errTable_te_s - errTable_te_sn
ptr - ptrn
pte - pten
ix_eff_s
ix_eff_sn

ax - axn



%--------------------------------
% Plot data (First 2 dimension)
%--------------------------------
% if data == 1,
%     figure,
%     subplot(2,2,1)
%     slr_view_data_multi(tte, xte, [1 2], ww_o)
%     axis equal;
%     title('SLR-Laplce version');
%     subplot(2,2,2)
%     slr_view_data_multi(tte, xte, [1 2], ww_f)
%     axis equal;
%     title('SLR-Variational version');
%     subplot(2,2,3)
%     slr_view_data_multi(tte, xte, [1 2], ww_s)
%     axis equal;
%     title('SMLR');
%     subplot(2,2,4)
%     slr_view_data_multi(tte, xte, [1 2], ww_r)
%     axis equal;
%     title('RMLR');
% end
% 
% fprintf('Finish demo !\n');
