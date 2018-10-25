% test iSLR 
%
% 2015/03/17
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.

clear
close all

data = 1;  % data = 1 or 2 or 3 (optional when you need to download testdata)

fprintf('This code demonstrates how iSLR works  ...\n');
%----------------------------
% Generate Data
%----------------------------
switch data
     case 1,
        D = 52;
        Ntr = 300;
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
        fprintf('Input feature dimension is %d. \n', D); 
        fprintf('The number of training samples is %d. \n', Ntr);

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
        fprintf('The mean value of the first 50 dimension is slightly different between two classes,\n');
        fprintf('whereas the remaining dimension has the same mean value.\n');
        fprintf('The degree of difference in the first 50 dimensions are manipulated\n')
        fprintf('by gradually changing the mean values in class 1 from 0 to 1.\n')
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
        fprintf('The data has already been preprocessed for classification.\n');
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
Algorithm = {'  SLR-VAR  ',...
             '  iSLR-VAR2 ',...
             '  iSLR-VAR3 '};

         
         
         
         
         for ii = 1 : length(Algorithm)
             tic
             fprintf('\n%s!!\n', Algorithm{ii});
             switch ii

                 case 1,


                     [ww_f, ix_eff_f, errTable_tr(:,:,ii), errTable_te(:,:,ii)] = biclsfy_slrvar(xtr, ttr, xte, tte,...
                         'nlearn', 300, 'mean_mode', 'each', 'scale_mode', 'each', 'invhessian',0);

                 case 2,

                     [ww_f2, ix_eff_f2, errTable_tr(:,:,ii), errTable_te(:,:,ii)] = biclsfy_islrvar(xtr, ttr, xte, tte, ii, ...
                         'nlearn', 300, 'mean_mode', 'each', 'scale_mode', 'each', 'invhessian',0);

                 case 3,

                     [ww_f3, ix_eff_f3, errTable_tr(:,:,ii), errTable_te(:,:,ii)] = biclsfy_islrvar(xtr, ttr, xte, tte, ii, ...
                         'nlearn', 300, 'mean_mode', 'each', 'scale_mode', 'each', 'invhessian',0);


             end


             time(ii)=toc;

         end





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
    

%-------------------------------
% Weight Parameters
%-------------------------------

figure, 
plot([ww_f],'b-','linewidth',3);
hold on
plot([ww_f2],'r-','linewidth',2);
plot([ww_f3],'g-','linewidth',1);
xlabel('Feature indicies');
ylabel('Weight')
legend(Algorithm)

%--------------------------------
% Plot data (First 2 dimension)
%--------------------------------

figure,
subplot(3,1,1)
slr_view_data(tte, xte, [1 2], ww_f(:,1))
axis equal;
title(Algorithm{1});
subplot(3,1,2)
slr_view_data(tte, xte, [1 2], ww_f2(:,1))
axis equal;
title(Algorithm{2});
subplot(3,1,3)
slr_view_data(tte, xte, [1 2], ww_f3(:,1))
axis equal;
title(Algorithm{3});




fprintf('\nFinish demo !\n');
