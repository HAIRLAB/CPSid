clear
clc
close all

addpath('./tools');
addpath('./data');
addpath('./SLR_dev');

%% design
np = 0;
txt_file = ['hc_n',num2str(np),'.txt'];
fileID = fopen(txt_file);



formatSpec = '%s';
N = 3;
C_text = textscan(fileID,formatSpec,N,'Delimiter','|');
C_data0 = textscan(fileID,'%f %f %f ');

num = 2000;

y = C_data0{1}(1:num,:);
u = C_data0{2}(1:num,:);
state = C_data0{3}(1:num,:);


%% make dic




A= [ones(size(u)) u u.^2  1./exp(u) u.^3./exp(u)  cos(2*u)./sin(u).^3];




%% identify subsystem
switch np
    case 0
        parameter.lambda = [0.005 0.005];   % the lambda of z in algorithm 1.
        parameter.epsilon = [1e-4  0.04];
    case 2
        parameter.lambda = [0.005 0.001];   % the lambda of z in algorithm 1.
        parameter.epsilon = [1e-4  0.04];
    case 4
        parameter.lambda = [0.008 0.005];   % the lambda of z in algorithm 1.
        parameter.epsilon = [1e-4  0.08];
    case 6
        parameter.lambda = [0.03 0.008];   % the lambda of z in algorithm 1.
        parameter.epsilon = [1e-4  0.105];
        
end

parameter.MAXITER = 5;
parameter.max_s = 20;%the max s
parameter.Phi = A;
parameter.y = y;
parameter.normalize_y = 1;
[result]=ihyde(parameter);

result.epsilon = parameter.epsilon(2);
result.lambda = parameter.lambda(2);
result.threshold = [0.05];
final_result  = finetuning( result);
sys = final_result.sys;
idx_sys = final_result.idx;

%%
ture_sys1 = [-0.5 ;1 ;0.5;0;0;0];
ture_sys2 = [0.5; 1 ;-0.5;0;0;0];

threshold1 = 0.05;%for no n2 n4
judge = ones(size(state));
final_idx{1} = [];
final_idx{2} = [];
final_idx{3} = [];
for k = 1:size(sys,2)
    
    if abs(sys(:,k)-ture_sys1)<threshold1 
        judge(idx_sys{k}) = 2;
    
    elseif abs(sys(:,k)-ture_sys2)<threshold1 
        judge(idx_sys{k}) = 3;

    end
    
end


wrong_position = find(judge - state ~=0 );

%% figure
% close all
% figure
% hold on
% plot(x,y,'r.','MarkerSize',20)
% xlabel('u','FontSize',13)
% ylabel('y','FontSize',13)
% 
% 
% ansy1 = zeros(size(y(:,1)));
% ansy1(find(judge==2)) = A(find(judge==2),:)*sys(:,1);
% ansy1(ansy1==0)=nan;
% plot(x,ansy1,'b.','MarkerSize',20)
% 
% ansy2 = zeros(size(y(:,1)));
% ansy2(find(judge==3)) = A(find(judge==3),:)*sys(:,2);
% ansy2(ansy2==0)=nan;
% plot(x,ansy2,'g.','MarkerSize',20)
% 
% 
% 
% legend('true y','Subysystem_1','Subysystem_2')
%%


Phi2 = [ones(size(u)) u  1./u.^2  cos(u)./sin(u).^3  ];

para_log.idx_sys = idx_sys;
para_log.beta = 10;

para_log.y = y;
para_log.Phi2 = Phi2;

para_log.normalize = 1;

[syslogic,labelMat,data] = ihydelogic(para_log);



