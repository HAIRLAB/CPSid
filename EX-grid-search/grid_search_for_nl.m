clear
clc
close all

addpath('../tools');
addpath('../data');

%% design


fileID = fopen('nl_n6.txt');

formatSpec = '%s';
N = 4;
C_text = textscan(fileID,formatSpec,N,'Delimiter','|');
C_data0 = textscan(fileID,'%f %f %f %f');
num = 4000;


y = C_data0{1}(1:num,:);
x1 = C_data0{2}(1:num,:);
x2 = C_data0{3}(1:num,:);

state = C_data0{4}(1:num,:);


%% make dic





A = [(x1+x2)./(x1-x2) x1./(6+x2)   x1.*x2  x1 x2 sin(x1) sin(x2) x1.^2 x2.^2];
para_num = 1;


train_index = 1:500;
y_train = y(train_index);
A_train = A(train_index,:);
check ={};
save_name  = 'cv_nl-all-jorge-1';

%% identify subsystem
for lam_z_i = 1:7
    lam_z = 1e-7*10^(lam_z_i);
    for lam_w_i = 1:3
        lam_w = 1e-4*10^(lam_w_i);
        for e_z_i = 1:2
            e_z = 1e-5*10^(e_z_i);
            for e_w_i = 1:2:9
                e_w = 1e-3*(e_w_i)*norm(y_train);
                par(para_num,:) = [lam_z lam_w e_z e_w];
                [lam_z lam_w e_z e_w]
                parameter.lambda = [lam_z lam_w];   % the lambda of z in algorithm 1.
                parameter.MAXITER = 5;
                parameter.max_s = 20;%the max s
                parameter.epsilon = [e_z e_w];
                parameter.Phi = A_train;
                parameter.y = y_train;
                parameter.normalize_y = 1;
                try
                [result]=ihyde(parameter);         
                result.epsilon = parameter.epsilon(2);
                result.lambda = parameter.lambda(2);
                result.threshold = [0.05];
                catch
                    result.early_stop =1;
                 
                end
                if result.early_stop == 1
                     check_tmp.fail = 1;
                     check{para_num} = check_tmp;
                     para_num = para_num+1;
                    continue;
                end
                final_result  = finetuning( result);
                sys = final_result.sys;
                idx_sys = final_result.idx;
                check_tmp.sys = sys;
                check_tmp.fail = 0;
                check_tmp.idx_sys = idx_sys;
                check_tmp.result = final_result;
                check_tmp.par = [lam_z lam_w e_z e_w];
                check{para_num} = check_tmp;
                para_num = para_num+1;
                
                
            end
            save(save_name);
        end
       save(save_name);
    end
     save(save_name);
end
save(save_name);
