%%
clear all
clc
close all

addpath('./tools');
addpath('./SLR_dev');
addpath('./data');

%%  load Data
load('ex_SmartGrid_model.mat')

% The nodes occur topology switching
choose_node=[11 12 22];
final_result = cell(0);

%% identify subsystem
for jj =1:size(choose_node,2)
    %% S is power data at the node, A is dictionary matrix
    S=[]; A=[];
    cur_node = choose_node(jj);
    for k=1:size(V,1)
        S = [S;-P(k,cur_node-1);-Q(k,cur_node-1)];
        A = [A;real(V(k,cur_node).*conj(V(k,:))) imag(V(k,cur_node).*conj(V(k,:)));...
            imag(V(k,cur_node).*conj(V(k,:))) -real(V(k,cur_node).*conj(V(k,:)))];
    end
    
    parameter.lambda = [5e-3 1e-6];  % the lambda of z in algorithm 1.
    parameter.MAXITER = 5;
    parameter.max_s = 20;%the max s
    parameter.epsilon = [0.015 0.05];
    parameter.Phi = A;
    parameter.y = S;
    parameter.normalize_y = 1;
    [result]=ihyde(parameter);
    
    result.epsilon = parameter.epsilon(2);
    result.lambda = parameter.lambda(2);
    result.threshold = [0.01];
    final_res  = finetuning( result);
    sys = final_res.sys;
    idx_sys = final_res.idx;
    final_result{jj} = final_res;
    
    %% identify logic
    index = 1:size(A,1);
    index = index';
    % the voltage diviation
    deltaV = [];
    deltaV_mag = [];
    for i=2:size(V,1)-1
        deltaV(i+1,:) = (V(i,:))-(V(i-1,:));
        deltaV_mag(i+1,:) = abs(V(i,:))-abs(V(i-1,:));
    end
    
    vol=[];
    for kk=1:33
        vol = [vol reshape([deltaV_mag(:,kk) deltaV_mag(:,kk)]', length(index),1)];
    end
    
    Phi2 = [vol ones(length(index),1) ];
    para_log.idx_sys = idx_sys;
    
    para_log.beta = 5.8e-5;
    para_log.y = S;
    para_log.Phi2 = Phi2;
    para_log.normalize = 1;
    [syslogic{jj},labelMat,data] = ihydelogic(para_log);
    
     T21(:,jj) = syslogic{1,jj}{2,1}/max(abs(syslogic{1,jj}{2,1}));
     T12(:,jj) = syslogic{1,jj}{1,2}/max(abs(syslogic{1,jj}{1,2}));

end
fprintf('The identified logic T12 is:')
display(T12)
fprintf('The identified logic T21 is:')
display(T21)