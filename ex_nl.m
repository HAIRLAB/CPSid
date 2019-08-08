clear
clc
close all

addpath('./tools');
addpath('./data');
addpath('./SLR_dev');
%% design

np = 0;
txt_file = ['nl_n',num2str(np),'.txt'];
fileID = fopen(txt_file);

formatSpec = '%s';
N = 4;
C_text = textscan(fileID,formatSpec,N,'Delimiter','|');
C_data0 = textscan(fileID,'%f %f %f %f');
num = 2000;


y = C_data0{1}(1:num,:);
x1 = C_data0{2}(1:num,:);
x2 = C_data0{3}(1:num,:);

state = C_data0{4}(1:num,:);


%% make dic

% A = [(x1+x2)./(x1-x2) 6*x1./(6+x2)   x1.*x2  ];

A = [(x1+x2)./(x1-x2) 6*x1./(6+x2)   x1.*x2  x1 x2 sin(x1) sin(x2) x1.^2 x2.^2];

%%
switch np
    case 0
        parameter.lambda = [1e-5 0.03];   % the lambda of z in algorithm 1.
        parameter.epsilon = [1e-4  0.6];
    case 2
        parameter.lambda = [5e-5 0.03];   % the lambda of z in algorithm 1.
        parameter.epsilon = [1e-4  0.8];
    case 4
        parameter.lambda = [1.5e-4 0.032];   % the lambda of z in algorithm 1.
        parameter.epsilon = [1e-4  2];
    case 6
        parameter.lambda = [1.5e-4 0.0282];   % the lambda of z in algorithm 1.
        parameter.epsilon = [1e-4  2];
        
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

threshold1 = 0.05;%for 0
threshold2 = 0.05;%for 1
judge = ones(size(state));
final_idx{1} = [];
final_idx{2} = [];
final_idx{3} = [];
for k = 1:size(sys,2)
    
    if abs(sys(1,k)-1)<threshold2 && abs(sys(2,k))<threshold2 && abs(sys(3,k))<threshold2
        judge(idx_sys{k}) = 2;
        
        
    elseif abs(sys(1,k))<threshold2 && abs(sys(2,k)-1)<threshold1 && abs(sys(3,k))<threshold2
        judge(idx_sys{k}) = 3;
        
    elseif abs(sys(1,k))<threshold2 && abs(sys(2,k))<threshold2 && abs(sys(3,k)-1)<threshold1
        judge(idx_sys{k}) = 4;
        
    end
    
end


wrong_position = find(judge - state ~=0 );

%% figure
% close all
% figure
% hold on
% plot(y,'r.','MarkerSize',20)
% xlabel('u','FontSize',13)
% ylabel('y','FontSize',13)
%
%
% ansy = zeros(size(y(:,1)));
% for k =1:size(sys,2)
%
%     ansy(idx_sys{k}) = A(idx_sys{k},:)*sys(:,k);
%
%
%
% end
%
% plot(ansy,'b','MarkerSize',20);


%
% ansy1 = zeros(size(y(:,1)));
% ansy1(find(judge==2)) = A(find(judge==2),:)*sys{1};
% ansy1(ansy1==0)=nan;
% plot(ansy1,'b.','MarkerSize',20)
%
% ansy2 = zeros(size(y(:,1)));
% ansy2(find(judge==3)) = A(find(judge==3),:)*sys{2};
% ansy2(ansy2==0)=nan;
% plot(ansy2,'g.','MarkerSize',20)
%
%
%
% legend('true y','Subysystem_1','Subysystem_2')


%%

Phi2 = [ones(size(x1)) x1 x2 exp(x1+x2)  x1.*x2 x1.^2 x2.^2];

para_log.idx_sys = idx_sys;
para_log.beta = 0.01;

para_log.y = y;
para_log.Phi2 = Phi2;

para_log.normalize = 1;
[syslogic,labelMat,data] = ihydelogic(para_log);

