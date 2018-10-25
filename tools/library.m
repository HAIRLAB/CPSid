function [yout ]= library(yin,polyorder,memory,basis_function)

nVars=size(yin,2);
cba = [];
if(memory>=1)
    n_time = size(yin,1);
    pretheta=[];
    
    for k =1:nVars
        abc=zeros(n_time,n_time-1);
       
        
        for i =1:n_time-1 
            n1=1;
            for j=1+i:n_time
                abc(j,i)=yin(n1,k);
                n1=n1+1;
            end
        end
        pretheta= [pretheta,abc];
        
    end
    
    
    pretheta=pretheta';
    
    cba=[];
    for i =1:n_time-1:(nVars)*(n_time-1)
        cba=[cba;pretheta(i:i+(memory-1),:)];
        
    end
    cba=cba';
    

end

yin = [ yin , cba];




nVars = size(yin,2);
n = size(yin,1);
% yout = zeros(n,1+nVars+(nVars*(nVars+1)/2)+(nVars*(nVars+1)*(nVars+2)/(2*3))+11);

ind = 1;
% poly order 0
yout(:,ind) = ones(n,1);
ind = ind+1;


% poly order 1
for i=1:nVars
    yout(:,ind) = yin(:,i);
    ind = ind+1;
end

if(polyorder>=2)
    % poly order 2
    for i=1:nVars
        for j=i:nVars
            yout(:,ind) = yin(:,i).*yin(:,j);
            ind = ind+1;
        end
    end
end

if(polyorder>=3)
    % poly order 3
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k);
                ind = ind+1;
            end
        end
    end
end

if(polyorder>=4)
    % poly order 4
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k).*yin(:,l);
                    ind = ind+1;
                end
            end
        end
    end
end

if(polyorder>=5)
    % poly order 5
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    for m=l:nVars
                        yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k).*yin(:,l).*yin(:,m);
                        ind = ind+1;
                    end
                end
            end
        end
    end
end


if  strcmp(basis_function.work,'on')
    nba = [];
    for i = basis_function.choose
        bf = [];
        for j = 1:size(yin,2)
            func =[];
            switch i
                case 1
                    for k = 1:size(basis_function.parameter_func1)
                        
                        func_k = basis_function.func1(yin(:,j),basis_function.parameter_func1(k,:));
                        func =[func func_k];
                    end
                case 2
                    for k = 1:size(basis_function.parameter_func2)
                        
                        func_k = basis_function.func2(yin(:,j),basis_function.parameter_func2(k,:));
                        func =[func func_k];
                    end
                case 3
                    for k = 1:size(basis_function.parameter_func3)
                        
                        func_k = basis_function.func3(yin(:,j),basis_function.parameter_func3(k,:));
                        func =[func func_k];
                    end
                    
                case 4
                    for k = 1:size(basis_function.parameter_func4)
                        
                        func_k = basis_function.func4(yin(:,j),basis_function.parameter_func4(k,:));
                        func =[func func_k];
                    end
                case 5
                    for k = 1:size(basis_function.parameter_func5)
                        
                        func_k = basis_function.func5(yin(:,j),basis_function.parameter_func5(k,:));
                        func =[func func_k];
                    end
            end
            bf =[bf,func];          
        end
        nba = [nba ,bf];
    end   
    yout = [yout , nba]; 
    
    
end




