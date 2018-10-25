function [ sys] = solver( parameter)

lambda = parameter.lambda;
A = parameter.A;
y = parameter.y;
y_normalize = parameter.y_normalize;
MAXITER = parameter.MAXITER;
if_print = parameter.if_print;
whole_index = parameter.whole_index;
if whole_index == 1
    idx_sys{1} = 1:length(y);
else
    idx_sys = parameter.idx_sys;
end
for k=1:size(idx_sys,2)
    if size(idx_sys{k},2) ==0
        continue;
    end
    T = A(idx_sys{k},:);
    for i = 1 : size(T,2)
        t(i) = norm(T(:,i),2);
        if t(i)==0
            t(i)=1;
        end
        T(:,i) = T(:,i)/t(i);
    end
    y1 = y(idx_sys{k});
    if y_normalize == 1;
        
        ty = max(abs(y1));
        if ty==0
            ty=1;
        end
        y1 = y1/ty;
    else ty = 1;
    end
        
    w = sparsesolver(y1, T,lambda , MAXITER,  parameter.Omega,if_print);
    sys(:,k) = w(:,end)./t'*ty;
end



end

