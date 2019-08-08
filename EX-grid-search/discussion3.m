clc 
clear

addpath('../tools');
basis_function.work='off';
%% success

dt = .005;
tspan=[dt:dt:0.5];
x0 = 0.99;
xall =[];
state_all =[];
for k=1:20
    
    N = length(tspan);
    if floor(k/2)==k/2
        [t,x] = ode45(@(t,x) -cos(x),tspan,x0);
        state = 2*ones(size(x));
    else
        
        [t,x] = ode45(@(t,x) - x^3,tspan,x0);
        state = 1*ones(size(x));
    end
    
    xall = [xall; x];
    state_all = [state_all; state];
    x0 = x(end);
    
end



x = xall;


% plot(x)
%%
state = state_all;
polyorder = 5;
memory = 0;
t = dt*1:length(x);
[x, dx, index] = estimatediff(x, t, 'solver', 1, []);
dx = dx/dt;
A= library(x,polyorder,memory,basis_function);
check ={};
save_name  = 'grid_search_d3-210';
index = 1:1000;
para_num = 1;
A_train = A(index,:);
y_train = dx(index);
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

