function result = Method1(parameter)
Thre = parameter.epsilon;
Phi = parameter.Phi;
[T,~] = size(Phi);
idx = (1:T)';
deleted_idx=[];

% secion 3.4 algorithm
for k = 1 : parameter.max_s
    if parameter.MAXITER>1
        fprintf(2,'the iteration is wrong!\n')
    else
        theta = ref18(parameter);
    end
    err=abs(parameter.y-parameter.Phi*theta)./(norm([1; theta],2)*sqrt(sum(([parameter.y -parameter.Phi]).^2,2)));
%     plot(err)
    delete_idx = find(err<Thre);
        
    if(k > 1)
        idx_sys{k} = setdiff(delete_idx, deleted_idx); %   only save the kth delete_idx
    else
        idx_sys{k} = delete_idx;
    end
%     idx = setdiff(idx,delete_idx);
    parameter.Phi(delete_idx,:) = 0;
    parameter.y(delete_idx,:) = 0;
    result.sys{k} = theta; 
    deleted_idx = find(parameter.y==0);
    if length(setdiff(idx,deleted_idx))<=30
        result.s = k ;
        fprintf(2,['There are ',num2str(k),' subsystems in total.\n'])
        break;
    end
end
result.idx = idx_sys;

end