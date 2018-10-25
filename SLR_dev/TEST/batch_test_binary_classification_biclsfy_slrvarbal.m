% 2 class classfication by linear SLR
% Simple classification when multiple features are relevant but either of
% them is not so strongly relevant.
% 
% discriminant fucntion : linear function 
% training data : simulated data generated from Gaussian Mixtures 
% test data     : simulated data generated from Gaussian Mixtures 

clear
close all

data = 2;
Nmonte = 50;
Nsenario = 2;

ERR_SLR = zeros(2,2,Nsenario);
ERR_BAL = zeros(2,2,Nsenario);


for ss = 1 : Nsenario

    switch ss
        case 1
            Ntrs = [90 10];
            Ntes = [1000 1000];
        case 2
            Ntrs = [90 10];
            Ntes = [1800 200];
        otherwise
            error('no such senario !');
    end

    for nn = 1 : Nmonte



        fprintf('This is a demo how the sparse logistic regression model works...\n');
        %----------------------------
        % Generate Data
        %----------------------------
        switch data
            case 1,

                D = 800;
                % mean
                mu1 = zeros(D,1);
                mu2 = [1.5; 0; zeros(D-2,1)];
                % covariance
                S = diag(ones(D,1));
                ro = 0.8;
                S(1,2) = ro;
                S(2,1) = ro;


                [ttr, xtr, tte, xte, g] = gen_simudata2([mu1 mu2], S, Ntrs, Ntes);

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
                D = 800;
                % mean
                mu1 = zeros(D,1);
                mu2 = [[1:-0.02:0]'; zeros(D-51,1)];
                % covariance
                S = diag(ones(D,1));

                [ttr, xtr, tte, xte, g] = gen_simudata2([mu1 mu2], S, Ntrs, Ntes);

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

        tic
        fprintf('\nNew Fast version (ARD-Variational)!!\n')
        [ww_fn, ix_eff_fn, errTable_tr_fn, errTable_te_fn, parm, AXall_fn, ptr_n, pte_n] = biclsfy_slrvar(xtr, ttr, xte, tte,...
            'nlearn', 300, 'mean_mode', 'none', 'scale_mode', 'none', 'invhessian', 0);
        toc

        PTE_SLR(nn,ss) = calc_percor(errTable_te_fn);
        ERR_SLR(:,:,ss) = ERR_SLR(:,:,ss) + errTable_te_fn;
        
        tic
        fprintf('\nNew Fast version (ARD-Variational)!!\n')
        [ww_fnc, ix_eff_fnc, errTable_tr_fnc, errTable_te_fnc, parm, AXall_fn, ptr_n, pte_n] = biclsfy_slrvarbal(xtr, ttr, xte, tte,...
            'nlearn', 300, 'mean_mode', 'none', 'scale_mode', 'none', 'balanced', 0, 'invhessian', 0);
        toc


        PTE_BAL(nn,ss) = calc_percor(errTable_te_fnc);
        ERR_BAL(:,:,ss) = ERR_BAL(:,:,ss) + errTable_te_fnc;
        
        
    end

end


for ss = 1 : Nsenario
   
    sensp_SLR = diag(ERR_SLR(:,:,ss))./sum(ERR_SLR(:,:,ss),2);
    sensp_BAL = diag(ERR_BAL(:,:,ss))./sum(ERR_BAL(:,:,ss),2);
    
    fprintf('--------- Senario %d  ---------------- \n', ss)
    fprintf('sensitivity  unbal : %0.3f  bal : %0.3f \n', sensp_SLR(1), sensp_BAL(1));   
    fprintf('specificity  unbal : %0.3f  bal : %0.3f \n', sensp_SLR(2), sensp_BAL(2));   
    
    
end



    


