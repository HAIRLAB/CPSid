clear
clc
close all

addpath('./tools');
addpath('./data');

%% design


load('current.mat')

%% make dic
state  = zeros(size(current));

state(1:19995) = 1;
t = 0:0.001:length(current)*0.001;

index = 1:300:length(current);
current = current(index);
state = state(index);
t1= t(index)';

y = current(2:end);
state = state(2:end);
t1= t1(2:end);
x = current(1:end-1);
polyorder = 2;

memory =5;
basis_function.work = 'off';
dic= library(x,polyorder,memory,basis_function);
% dic(:,2) = [];
A = dic(memory+1:end,:);
y = y(memory+1:end,:);
t1 = t1(memory+1:end,:);
state = state(memory+1:end,:);
A111 = A;
y111 = y;

%% identify subsystem

% 
% parameter.lambda = [1.5 0.00001];   % 65
% parameter.MAXITER = 5;
% parameter.max_s = 20;%the max s
% parameter.epsilon = [1e-4  0.026];


parameter.lambda = [1.5 0.00005];   % 64
parameter.MAXITER = 5;
parameter.max_s = 20;%the max s
parameter.epsilon = [1e-4  0.026];

parameter.Phi = A;
parameter.y = y;
parameter.normalize_y = 1;
[result]=ihyde(parameter);





result.epsilon = parameter.epsilon(2);
result.lambda = parameter.lambda(2);
result.threshold = [0.05];

result.mode = 'knn';
result.half_win_size = 5;
final_result  = finetuning( result);
sys = final_result.sys;
idx_sys = final_result.idx;



%%


% % clear syslogic labelMat data para_log Phi2
% Phi2 = [ones(size(x)) x ];
t = 1:length(y);

Phi2 = [ones(size(t')) t' ];

para_log.idx_sys = idx_sys;
para_log.beta = 0.1;

para_log.y = y;
para_log.Phi2 = Phi2;

para_log.normalize = 1;

[syslogic,labelMat,data] = ihydelogic(para_log);


%%
close all
figure
hold on
index1 = 1:length(state);
y11 = y(index1);
judge11 = state(index1);
x11 = x(index1);
A1 = A(index1,:);

ansy1 = 1000000*ones(size(y11(:,1)));
ansy1(find(judge11==0)) = y11(find(judge11==0));
ansy1(ansy1==1000000)=nan;
plot(ansy1,'b','Linewidth',5);

ansy11 = 1000000*ones(size(y11(:,1)));
ansy11(idx_sys{2},:) = A1(idx_sys{2},:)*sys(:,2);
ansy11(ansy11==1000000)=nan;
plot(ansy11,'bo','MarkerSize',10)
err1 = ansy1 - ansy11;

ansy2 = 1000000*ones(size(y11(:,1)));
ansy2(find(judge11==1)) = y11(find(judge11==1));
ansy2(ansy2==1000000)=nan;
plot(ansy2,'r','Linewidth',5)

ansy22 = 1000000*ones(size(y11(:,1)));
ansy22(idx_sys{1},:) = A1(idx_sys{1},:)*sys(:,1);
ansy22(ansy22==1000000)=nan;
plot(ansy22,'ro','MarkerSize',10)
err2 = ansy2 - ansy22;
