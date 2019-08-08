function [result]=ihyde(parameter)

%
% -------------------------------------------------------------------------
% INPUT parameter
% -------------------------------------------------------------------------
% parameter.Phi: constructed dictionary matrix phi. Each row is a datapoint.
%
% parameter.y: column vector containing the output datapoints.
%
% parameter.max_s: the max number of subsystems that could be identified by IHYDE
%
% parameter.lambda: a 1 by 2 vector(lambda_z, lambda_w). For the identification of a new subsystems.
%
% parameter.epsilon: a 1 by 2 vector(epsilon_z, epsilon_w). For the identification of a new subsystems.
%
% parameter.MAXITER: the max number of iterations that the sparsesolver function solves.
%
% parameter.normalize_y : set to 1 if y need to be normalized.
%
%
% -------------------------------------------------------------------------
% OUTPUT result
% -------------------------------------------------------------------------
%
% result.sys : the model of each subsystem.
%
% result.idx_sys :the index of each subsystem.
%
% result.theta: z of each iter.
%
% result.error : each subsystem's fitting error.
%
% result.deleted_idx: index of the already deleted data points
%
% result.no_delete_idx: the rest data points
%

%
% Copyright is with the following author:
%
% (C) 2016 Ye Yuan,
Phi = parameter.Phi;
lambda = parameter.lambda;
epsilon = parameter.epsilon;
MAXITER = parameter.MAXITER;
y = parameter.y;
result = parameter;
deleted_idx = [];
result.early_stop = 0;
for k = 1 : parameter.max_s

    T = null(Phi')';

    y1 = T*y;

    Omega = (T*T')^(-1);

    solver_parameter.y = y1;
    solver_parameter.A = T;
    solver_parameter.lambda = lambda(1);
    solver_parameter.MAXITER = MAXITER;
    solver_parameter.Omega = Omega;
    solver_parameter.y_normalize = 1;
    solver_parameter.if_print = 1;
    solver_parameter.whole_index = 1;
    fprintf(2,'Finding %dth subsystem\n', k);
    z = solver( solver_parameter);
    n =  size(y,1);

    delete_idx = [];

    result.theta(:,k) = abs(z);    %save theta in check for tunning lambda and epsilon
    for i = 1 : n
        if abs(z(i)) < epsilon(1)  %% compared with error and choose theta
            delete_idx = [delete_idx i];
        end
    end
    
    
    %
    if(k > 1)
        idx_sys{k} = setdiff(delete_idx, deleted_idx'); %   only save the kth delete_idx
        
    else
        idx_sys{k} = delete_idx;
    end
    
    
    if size(idx_sys{k},2) == 0
        result.early_stop = 1;
        result.s = k;
        break;
    end
    [delete_idx,    alone]= find_idx( idx_sys{k},0);%Find the longest continuous sequence
    result.longgest_idx{k} = delete_idx;
    
    T = Phi(delete_idx,:);
    solver_parameter.y = y(delete_idx);
    solver_parameter.A = T;
    solver_parameter.lambda = lambda(2);
    solver_parameter.MAXITER = MAXITER;
    solver_parameter.y_normalize = 1;
    solver_parameter.Omega = 1;
    solver_parameter.if_print = 1;
    solver_parameter.whole_index = 1;
    fprintf(2, 'Identifing the %dth subsystem\n',k);
    sys(:,k) = solver( solver_parameter);

    result.sys(:,k) = sys(:,k);
    error = abs(y - Phi * sys(:,k));
    result.sys_allerror(:,k) = abs(parameter.y - parameter.Phi * sys(:,k));
    result.error(:,k) = error;
    delete_idx = [];
    for i = 1 : n
        if error(i) < epsilon(2)  %%compared with error,it is more convenient to find a better epsilon by choosing  theta
            delete_idx = [delete_idx i];
        end
    end
    
    %
    if(k > 1)
        idx_sys{k} = setdiff(delete_idx, deleted_idx'); %   only save the kth delete_idx
        
    else
        idx_sys{k} = delete_idx;
    end
    
    result.idx(k,1:size(idx_sys{k},2)) = idx_sys{k}; %save kth idx_sys
    result.idx_sys{k} = idx_sys{k};
    
    y(delete_idx) = 0;
    Phi(delete_idx, :) = 0;
    
    
    no_delete_idx = find(y~=0); %save all the deleted_idx
    result.no_delete_idx = no_delete_idx;
    
    if size(no_delete_idx',2)~=0
        [~,    alone]= find_idx( no_delete_idx',0);%if alone = 1,it mearns the length of the longest consecutive sequence is less than the threshold
    end
    
    result.deleted_idx = find(y==0);
    
    deleted_idx = find(y==0); %save all the deleted_idx
    size(deleted_idx )
    result.deleted_size(k) = size(deleted_idx ,1);
    fprintf( 2,'After the identification of %dth subsystem,there are %d elements have been deleted.\n\n\n\n',k,size(deleted_idx,1 ));
    if size(deleted_idx,1 )==size(y,1) || (k>1&&result.deleted_size(k-1) == size(deleted_idx ,1))||alone==1%if all idx have been deleted ,break the loop and save the number of subsystems
        
        result.s = k;
        break;
    end

    result.s = k;
end



