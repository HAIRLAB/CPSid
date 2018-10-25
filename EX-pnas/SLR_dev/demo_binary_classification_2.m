% Binary classfication demos
% discriminant fucntion : linear function 
%
% updated 2010/03/15  
% * check whether optimization toolbox is avaliable 
% * 'scale_mode' and 'mean_mode' are changed. 
% updated 2009/06/15 
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.

clear
close all

data = 3;  % data = 1 or 2 or 3 (optional when you need to download testdata)

fprintf('This code demonstrates how a binary classification problem is solved ...\n');
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
        mu2 = [4; 0; zeros(D-2,1)];
        % covariance
        S = diag(ones(D,1));
        ro = 0.5;
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


Algorithm = {'  SLR-LAP  ',...
    '  SLR-VAR  ',...
    '    RVM    ',...
    '  RLR-VAR  ',...
    'L1-SLR-LAP ',...
    'L1-SLR-COMP'};


tic
fprintf('\n%s!!\n', Algorithm{2});
[ww_f, ix_eff_f, errTable_tr(:,:,2), errTable_te(:,:,2), parm, AXall,Ptr,Pte, XI2, EDF] = biclsfy_slrvar(xtr, ttr, xte, tte,...
    'nlearn', 100, 'nstep', 1, 'mean_mode', 'none', 'scale_mode', 'none', 'invhessian',1);
time(2)=toc;
toc


%
%  plot
%
figure
for nlearn = 1 : size(XI2,2),
    
    for jj = 1 : size(XI2,1)
        if ttr(jj) == 0
            plot(xtr(jj,1), xtr(jj,2), 'o', 'markersize', 1/XI2(jj, nlearn)*10);
        else
            plot(xtr(jj,1), xtr(jj,2), 'ro', 'markersize', 1/XI2(jj, nlearn)*10);
        end
        hold on
        
    end
    pause(0.1);
    if nlearn ~= size(XI2,2) clf; end
end
        
    

nx = [1:30];

figure,
for ii = 1 : 30
subplot(6,5,ii)
[ax,h1,h2]=plotyy(nx,log10(AXall(ii,nx+1)), nx, EDF(ii,nx));    
axis(ax(2), [nx(1) nx(end), -10 2]);
axis(ax(1), [nx(1) nx(end), -inf inf]);
set(h2, 'marker', '.');
set(h1, 'marker', '.');

end

figure,
for ii = 1 : 30
subplot(6,5,ii)
[ax,h1,h2]=plotyy(nx,log10(AXall(ii+30,nx+1)), nx, EDF(ii+30,nx));    
axis(ax(2), [nx(1) nx(end), -10 2]);
axis(ax(1), [nx(1) nx(end), -inf inf]);
set(h2, 'marker', '.');
set(h1, 'marker', '.');

end

figure
plot(sum(EDF,1))



