clear
clc
close all

addpath('./tools');
addpath('./data');
%%  load Data
data = load('TwoMissLine.mat');
data1 = data.allData{1};
data2 = data.allData{2};
u1 = data1{1}.';
i1 = data1{2}.';
y1 = data1{3}.';

u2 = data2{1}.';
i2 = data2{2}.';
y2 = data2{3}.';

index = 1:10;
u2 = u2(index,:);
i2 = i2(index,:);


%% library
A1=[];
for jj=1:size(u1,1)
    A1=[A1;real(u1(jj,:)) -imag(u1(jj,:)); imag(u1(jj,:)) real(u1(jj,:))];
end
A2=[];
for jj=1:size(u2,1)
    A2=[A2;real(u2(jj,:)) -imag(u2(jj,:)); imag(u2(jj,:)) real(u2(jj,:))];
end

% A1 = [real(u1) -imag(u1);imag(u1) real(u1) ];
% A2 = [real(u2) -imag(u2);imag(u2) real(u2) ];

A = [A1;A2];

%i = [ii1;ii2];
%% identify subsystem
choose = 6;
ii1 = reshape(([real(i1(:,choose)) imag(i1(:,choose))])',length(i1(:,1))*2,1);
ii2 = reshape(([real(i2(:,choose)) imag(i2(:,choose))])',length(i2(:,1))*2,1);
i = [ii1; ii2];


parameter.lambda = [1e-3 1e-6];  % the lambda of z in algorithm 1.
parameter.MAXITER = 5;
parameter.max_s = 20;%the max s
parameter.epsilon = [0.008  0.05];



parameter.Phi = A;
parameter.y = i;
parameter.normalize_y = 1;
[result]=ihyde(parameter);


result.epsilon = parameter.epsilon(2);
result.lambda = parameter.lambda(2);
result.threshold = [0.05];
final_result  = finetuning( result);
sys = final_result.sys;
idx_sys = final_result.idx;

%% identify logic
index = 1:size(A,1)/2;
index = kron(index,ones(1,2));
index = index';
Phi2 = [ones(length(index),1) index i i.^2 index./(index+10*i)];


para_log.idx_sys = idx_sys;
para_log.beta = 0.001;

para_log.y = i;
para_log.Phi2 = Phi2;

para_log.normalize = 1;
[syslogic,labelMat,data] = ihydelogic(para_log);
