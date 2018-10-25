function [data_new, derivative, index] = estimatediff(data, t, type, kth, order)

% Estimate the derivative
% version 1.0 : we only estimate the first derivative, i.e. order = 1.
%            we will add higher order estimation in the future version
% please cite Kris De Brabanter's JMLR paper: Derivative Estimation with Local Polynomial Fitting



%%
% clear all; close all
% 
% type = 'solver';
% %                 type = 'central';
% %                 type = 'eularforward';
% %                 type = 'eularbackward';
% tt= 0:pi/10:24;
% y1 = tt.*sin(tt/2)-tt;
% y1 = y1(:);
% y2 = cos(tt);
% y2 = y2(:);
% data = [y1, y2];
% noise =0*randn(size(data));
% data = data+noise;
% derivative_true1 = sin(tt/2) + tt.*cos(tt/2)/2-1;
% derivative_true2 = -sin(tt);
% derivative_true = [derivative_true1(:), derivative_true2(:)];
% t = tt;
% which_var  = 1;
% kth =2; % please specify the k value


%% ================== solver =======================
plot_onoff='off';

n = size(data,1);%n是data的行数

if ~exist('t') || isempty(t)%~表示将exist(t)的结果取反，isempty(t)当t为空矩阵时返回1
    t = 1:n;
    t = t(:);  %t是列向量，从1到n
end

if ~isempty(t)
    t = t(:); %t不为空矩阵时，给其赋值
    if length(t)~=n
        error('dimension of t is not consistent');%若赋值后t的长度不为n，报错
    end
end

if ~exist('type') || isempty(type)
    type = 'central';
end

if strcmp(type, 'central')  %strcmp函数比较两个字符串，前者大于后者时返回正数，相等时返回0
    type = 'solver';
    kth = 1;
end

if ~exist('order') || isempty(order)
    order = 1;
end


switch type
    case 'solver'
        %%
        index = kth+1:n-kth;%k+1<i<n-k，i取值范围
        
        for j  = 1:1:kth
            w(j) = (6*j^2)/(kth*(kth+1)*(2*kth+1));  %构造wj
        end
        
        t_new = t(index,:);
        
       % for i = 1: size(data_new,1)
        for i = index
            for j = 1:kth
                diff_i(j,:) = w(j) * ...
                    (data(i+j,:) - data(i-j,:))./(t(i+j,:) - t(i-j,:));  %对称差分商diff_i
            end
            derivative(i,:) = sum(diff_i,1);  %对矩阵每一列的元素求和,得到一个行向量，即derivative的第i行
            d_forward(i,:) = (data(i+1,:) - data(i,:))./(t(i+1,:) - t(i,:));
            d_backward(i,:) = data(i,:) - data(i-1,:)./(t(i,:) - t(i-1,:));
        end
        derivative(1:kth,:) = [];%消除第1到kth行的数据
        d_forward(1:kth,:) = [];
        d_backward(1:kth,:) = [];
        
    case 'eularforward'
        %%
        index  = 1: n-1;
        for i = index
            d_forward(i,:) = (data(i+1,:) - data(i,:))./(t(i+1,:) - t(i,:));
        end
        derivative = d_forward(index,:);
        
    case 'eularbackward'
        %%
        index  = 2: n;
        for i = index
            d_backward(i,:) = data(i,:) - data(i-1,:)./(t(i,:) - t(i-1,:));
        end
        derivative = d_backward(index,:);
        index  = index - 1;
end

index = index(:);  %横行变成竖行
data_new = data(index,:);  %data_new是data的第k+1到n-k行

if size(data_new,1) ~= size(derivative,1) %检查维数是否一致
    error('please check dimensions')
end


%% to comment
if strcmp(plot_onoff,'on')
    derivative_true = derivative_true(index,:);
    % derivative_true = normalize_var(derivative_true, -1, 1);
    % derivative = normalize_var(derivative, -1, 1);
    % d_forward = normalize_var(d_forward, -1, 1);
    
    
    figure; plot(data_new(:,which_var),'r');
    figure; plot(derivative(:,which_var),'b-*'); hold on; plot(derivative_true(:,which_var),'r-o');
    hold on; plot(d_forward(:,which_var),'g')
    legend('estimated derivative','true derivative','eular')
end