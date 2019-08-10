clear all,close all,clc

addpath('./tools');
addpath('./data');
addpath('./SLR_dev');
%%  load Data


basis_function.work='off';
data  = load('normal_car.mat');  %%  load Data

flag = data.flag;   % 1:straghtway   0:curve
du = data.dy;
v = data.v/10;


%% library
polyorder = 4;
memory = 4;

A= library(v,polyorder,memory,basis_function);
v_k1 = v(memory+1:end-1,:);
v_k2 = v(memory:end-2,:);
v_k3 = v(memory-1:end-3,:);
v_k4 = v(memory-2:end-4,:);
v = v(memory+2:end,:);

A = A(memory+2:end,:);
du = du(memory+2:end);
flag = flag(memory+2:end);
%A = A(:,1:3);%withoutredundant dictionary functions.

%% identify subsystem

parameter.lambda = [0.01 1e-5];   % the lambda of z in algorithm 1.
parameter.MAXITER = 5;
parameter.max_s = 20;%the max s
parameter.epsilon = [100  8];
parameter.Phi = A;
parameter.y = du;
parameter.normalize_y = 1;
[result]=ihyde(parameter);



result.epsilon = parameter.epsilon(2);
result.lambda = parameter.lambda(2);
result.threshold = [0.05];
final_result  = finetuning( result);
sys = final_result.sys;
idx_sys = final_result.idx;


%% identify logic
Phi2 = [ones(size(flag)) flag  sin(v) cos(v) tan(v) (v_k1-v_k4)./v_k2  v_k1.*tan(v_k3)];
% Phi2 = [ones(size(flag)) flag ];

para_log.idx_sys = idx_sys;
para_log.beta = .001;
para_log.y = du;
para_log.Phi2 = Phi2;
para_log.normalize = 1;

[syslogic,labelMat,data] = ihydelogic(para_log);

%% compare with the right answer
judge = 7*ones(size(flag));
judge(idx_sys{1}) = 0;
judge(idx_sys{2}) = 1;
wrong_position = find((judge-flag)~=0);
ans_sys_idx{1} = find(flag==1);
ans_sys_idx{2} = find(flag==0);
%%
close all
figure(1)
axes1 = axes('Parent',figure(1));
hold on
color = {'r' ,'b'};
for i =1:2
    input1 = zeros(size(v));
    input1(ans_sys_idx{i},1) = v(ans_sys_idx{i},1);
    
    input1(input1==0)=nan;
    
    plot(input1(:,1),'Color',color{i},'LineWidth',3);
    
    
end
legend('Subsystem_1','Subsystem_2')

xlabel('Time(4ms)','FontWeight','bold');
ylabel('Velocity(cm/s)','FontWeight','bold');
box(axes1,'on');
set(axes1,'FontSize',14,'FontWeight','bold','LineWidth',1.5);
legend(axes1,'show');
%%
figure(2)
axes1 = axes('Parent',figure(2));
hold on
color = {'r' ,'b'};
input1 = zeros(size(du,1),2);
for i =1:2
    
    input1(ans_sys_idx{i},i) = du(ans_sys_idx{i},1);
    
    input1(input1(:,i)==0,i)=nan;
    
    plot(input1(:,i),'Color',color{i},'LineWidth',3);
    
    
end



legend('Subsystem_1','Subsystem_2')

xlabel('Time(4ms)','FontWeight','bold');
ylabel('du','FontWeight','bold');
box(axes1,'on');
set(axes1,'FontSize',14,'FontWeight','bold','LineWidth',1.5);
legend(axes1,'show');
