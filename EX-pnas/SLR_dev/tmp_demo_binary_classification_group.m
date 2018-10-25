% Binary classfication demos
% discriminant fucntion : linear function 
%
% last updated 2009/06/15 
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.

clear
close all

data = 2;  % data = 1 or 2 or 3 (optional when you need to download testdata)

fprintf('This code demonstrates how a binary classification problem is solved ...\n');
%----------------------------
% Generate Data
%----------------------------
switch data
     case 1,
        D = 52;
        Ntr = 30000;
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
        fprintf('But only the first dimension has difference in mean value between two classes,\n');
        fprintf('and the other dimension has same mean value.\n');
        fprintf('Therefore only the first dimension is detected as a meaningful feature,\n')
        fprintf('if you select features by feature-wise t-value ranking method.\n');
        fprintf('However due to the correlataion between the second dimension and the first dimension,\n')
        fprintf('inclusion of the second dimension makes classfication more accurate.\n');
        fprintf('For comparison, this demo also computes the classification performance of linear RVM, \n');
        fprintf('which is Bayesian counterpart of support vector machine (SVM).\n');
        fprintf('Input feature dimension is %d. \n', D); 
        fprintf('The number of training samples is %d. \n', Ntr);

    case 2,
        D = 100;
        D_eff = 3;
        Ntr = 100;
        Nte = 100;
        % mean
        mu1 = zeros(D,1);
        mu2 = [[1:-1/D_eff:0]'; zeros(D-D_eff-1,1)];
        % covariance
        S = diag(ones(D,1));

        [ttr, xtr, tte, xte, g] = gen_simudata([mu1 mu2], S, Ntr, Nte);

        fprintf('\nThe data is generated from 2 Gaussian Mixture model of which centers (mean) are different.\n');
        fprintf('The mean value of the first 50 dimension is slightly different between two classes,\n');
        fprintf('whereas the remaining dimension has the same mean value.\n');
        fprintf('The degree of difference in the first 50 dimensions are manipulated\n')
        fprintf('by gradually changing the mean values in class 1 from 0 to 1.\n')
        fprintf('For comparison, this demo also computes the classification performance of linear RVM, \n');
        fprintf('which is Bayesian counterpart of support vector machine (SVM).\n');
        fprintf('Input feature dimension is %d. \n', D); 
        fprintf('The number of training samples is %d. \n', Ntr);
        
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
        fprintf('Input feature dimension is %d. \n', D); 
        fprintf('The number of training samples is %d. \n', Ntr);
end

%--------------------------------
% Plot data (First 2 dimension)
%--------------------------------
slr_view_data(ttr, xtr);
axis equal;
title('Training Data')

fprintf('\n\nPress any key to proceed \n\n');
pause

%--------------------------------
% Learn Paramters
%--------------------------------
Algorithm = {'  SLR-LAP-Group10',...
             '  SLR-LAP-Group20',...
             '  SLR-LAP        ',...
       };
% 

for ii = 1 : 10 
group10((ii-1)*10+1:ii*10) = randn(1);
end
group10 = group10(:);

for ii = 1 : 20 
group20((ii-1)*5+1:ii*5) = ii;
end
group20 = group20(:);


tic
fprintf('\n%s!!\n', Algorithm{1});
[ww_og10, ix_eff_og10, errTable_tr(:,:,1) errTable_te(:,:,1), parm, AXog10] = biclsfy_slrlap(xtr, ttr, xte, tte,...
    'wdisp_mode', 'off', 'nlearn', 500, 'mean_mode', 'none', 'scale_mode', 'none', 'group', group10);
time(1) = toc;
toc

tic
fprintf('\n%s!!\n', Algorithm{2});
[ww_og20, ix_eff_og20, errTable_tr(:,:,2) errTable_te(:,:,2), parm, AXog20] = biclsfy_slrlap(xtr, ttr, xte, tte,...
    'wdisp_mode', 'off', 'nlearn', 500, 'mean_mode', 'none', 'scale_mode', 'none', 'group', group20);
time(2) = toc;
toc

fprintf('\n%s!!\n', Algorithm{3});
[ww_o, ix_eff_o, errTable_tr(:,:,3) errTable_te(:,:,3), parm, AXo] = biclsfy_slrlap(xtr, ttr, xte, tte,...
    'wdisp_mode', 'off', 'nlearn', 500, 'mean_mode', 'none', 'scale_mode', 'none');
time(3) = toc;
toc





%--------------------------------
%  Result Table
%--------------------------------
fprintf('\n\n')
fmt1 = sprintf(' %%%ds       %%3.2f      %%3.2f      %%3.3f \n',11);
fmt2 = sprintf('                 %%5s    %%5s   %%6s\n');
fprintf(fmt2, 'Train(%)','Test(%)','Time(sec.)');
for ii = 1 : length(Algorithm),
    fprintf(fmt1,  Algorithm{ii}, calc_percor(errTable_tr(:,:,ii)), calc_percor(errTable_te(:,:,ii)), time(ii));
end
    



%--------------------------------
% Plot data (First 2 dimension)
%--------------------------------

figure,
subplot(3,1,1)
slr_view_data(tte, xte, [1 2], ww_og10(:,1))
axis equal;
title('SLR-Laplce group10');
subplot(3,1,2)
slr_view_data(tte, xte, [1 2], ww_og20(:,1))
axis equal;
title('SLR-Laplce group20');
subplot(3,1,3)
slr_view_data(tte, xte, [1 2], ww_o(:,1))
axis equal;
title('SLR-Laplace');

figure,
subplot(1,3,1)
imagesc(log10(AXog10), [0 8])
title('SLR-Laplce group10');
subplot(1,3,2)
imagesc(log10(AXog20), [0 8])
title('SLR-Laplce group20');
subplot(1,3,3)
imagesc(log10(AXo), [0 8])
title('SLR-Laplace');
set(gcf,'name', 'Learning curve of relevance parameters');


fprintf('\nFinish demo !\n');
