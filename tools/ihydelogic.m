function [syslogic,labelMat,data,num] = ihydelogic(para_log)

beta = para_log.beta ;

idx_sys = para_log.idx_sys ;
s = size(idx_sys,2);
y = para_log.y;
Phi2 = para_log.Phi2 ;
normalize = para_log.normalize ;

t=1;
if normalize>1
    for i = 1 : size(Phi2,2)
        t(i) = norm(Phi2(:,i),2);
        if t(i)==0
            t(i)=1;
        end
        Phi2(:,i) = Phi2(:,i)/t(i);
    end
end




for k = 1 : s
    for k1 = 1 : s
        if(k1~=k)
            
            n =  size(y,1);
            
            
            
            idx1 = intersect(idx_sys{k}, idx_sys{k}+1);
            idx2 = intersect(idx_sys{k1}, idx_sys{k}+1);
            
            
            
            idx = [idx1 idx2];
            
            set{k,k1} = zeros(n,1);
            set{k,k1}(idx2) = 1;
            num(k,k1) = length(idx2);
            
            
            data{k,k1} =Phi2(idx,:);
            
            labelMat{k,k1} = set{k,k1}(idx,:);
            
            
            if (labelMat{k,k1}==0)
                coe{k,k1} = [];
            else
                
                
                [coe{k,k1}, ix_eff, W, AX] = slr_learning_l1(labelMat{k,k1}, data{k,k1}, @linfun, beta);
                coe{k,k1} = coe{k,k1}./t';
                
                %
            end
            
        end
    end
end
syslogic = coe;
end
