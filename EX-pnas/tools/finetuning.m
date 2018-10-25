function [ final_result ] = finetuning( result)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fprintf(1,'Finetuning each subsystems ...\n');
sys = result.sys;
idx_sys = result.idx_sys;
y = result.y;
A = result.Phi;
threshold = result.threshold;

epsilon = result.epsilon;


sys_mark = zeros(size(sys,2),1);
num = 0;
for i = 1:size(sys,2)
    
    if sys_mark(i) == 1
        continue;
    end
    idx_tmp = [];
    for j = i:size(sys,2)
        
        judge = abs(sys(:,i)-sys(:,j))<threshold;
        if min(judge) == 1
            sys_mark(j) = 1;
            idx_tmp = [idx_tmp idx_sys{j}];
            
        end
    end
    num = num+1;
    final_result.idx_sys{num} = idx_tmp;
end


result.Omega = 1;
final_result.y = result.y;
final_result.A = result.Phi;
final_result.lambda = result.lambda;

final_result.MAXITER = result.MAXITER;
final_result.Omega =  result.Omega;

for iter = 1:20
    iter

    final_result.y_normalize = 1;

    final_result.if_print = 0;
    final_result.whole_index = 0;
    sys = solver( final_result);
    allerror = [];
    for i = 1:size(sys,2)
        allerror(:,i) = abs(y - A*sys(:,i));
    end
    
    error_square = 0;
    for i = 1:size(y,1)
        error_square = min(allerror(i,:)).^2 + error_square;
    end
    
    mse = error_square/size(y,1);
    
    if iter>1&&(last_mse - mse)<=0
        break;
    end
    
    last_mse = mse;
    final_result.sys = sys;
    
    minerror_pos = zeros(size(y));
    
    for i = 1:size(result.deleted_idx,1)
        tmp = allerror(result.deleted_idx(i),:);
        if min(tmp)<=epsilon
            minerror_pos(result.deleted_idx(i)) = find(tmp==min(tmp));
        end
    end
    
    for i = 1:size(sys,2)
        final_result.idx{i} = find(minerror_pos==i);
    end
    
end
allerror = [];
for i = 1:size(final_result.sys,2)
    allerror(:,i) = abs(y - A*final_result.sys(:,i));
end
minerror_pos = [];
for i =1:size(y,1)
    minerror_pos(i) = find(allerror(i,:)==min(allerror(i,:)));
    
end
final_result.allerror = allerror;
final_result.minerror_pos = minerror_pos;
for i = 1:size(sys,2)
    final_result.idx{i} = find(minerror_pos==i);
end




