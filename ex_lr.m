clear
clc
close all

addpath('./tools');
addpath('./data');
addpath('./SLR_dev');

%% design
np = 0;
txt_file = ['lr_n',num2str(np),'.txt'];
fileID = fopen(txt_file);

formatSpec = '%s';
N = 7;
C_text = textscan(fileID,formatSpec,N,'Delimiter','|');
C_data0 = textscan(fileID,'%f%f%f%f%f%f%f');



start = 1;
num = 2000;

y = C_data0{1}(start:num);
x1 = C_data0{2}(start:num);
x2 = C_data0{3}(start:num);
x4 = C_data0{4}(start:num);
x3 = C_data0{5}(start:num);
x5 = C_data0{6}(start:num);
state = C_data0{7}(start:num);

%% make dic

A = [ones(size(y)) x1-x2  , 1./(x1-x2),x1.^2,x2.^2];
% A = [ x1-x2  , 1./(x1-x2)];

%%
switch np
    case 0
        parameter.lambda = [1e-4 1e-3];   % the lambda of z in algorithm 1.
        parameter.epsilon = [1e-4  0.005];
    case 2
        parameter.lambda = [1e-4 0.1];   % the lambda of z in algorithm 1.
        parameter.epsilon = [1e-4  0.06];
    case 4
        parameter.lambda = [5e-4 0.15];   % the lambda of z in algorithm 1.
        parameter.epsilon = [1e-4  0.2];
    case 6
        parameter.lambda = [1e-3 0.15];   % the lambda of z in algorithm 1.
        parameter.epsilon = [1e-4  0.2];
        
end
parameter.MAXITER = 5;
parameter.max_s = 20;%the max s

parameter.Phi = A;
parameter.A_norm = 1;
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
threshold1 = 0.05;%for compare with 1
threshold2 = 0.05;%for compare with 0
judge = ones(size(state));
final_idx{1} = [];
final_idx{2} = [];
final_idx{3} = [];
for k = 1:size(sys,2)
    
    if abs(sys(1,k))<threshold2 && abs(sys(2,k)+1)<threshold1 && abs(sys(3,k))<threshold2
        judge(idx_sys{k}) = 2;
        
        
    elseif abs(sys(1,k))<threshold1 && abs(sys(2,k))<threshold2 && abs(sys(3,k)-1)<threshold1
        judge(idx_sys{k}) = 3;
        
        
    elseif abs(sys(1,k))<threshold2 && abs(sys(2,k))<threshold2 && abs(sys(3,k))<threshold2
        judge(idx_sys{k}) = 4;
        
    end
    
end


wrong_position = find(judge - state ~=0 );


%% figure
% close all
% figure(start)
% hold on
%
% plot(input,y,'r.','MarkerSize',20)
% xlabel('u','FontSize',13)
% ylabel('y','FontSize',13)
%
%
%
% if size(find(judge==2),1)>0
% ansy1 = zeros(size(y(:,1)));
% ansy1(find(judge==2)) = A(find(judge==2),:)*sys(:,1);
% ansy1(ansy1==0)=nan;
% plot(input,ansy1,'b.','MarkerSize',20)
% end
%
% if size(find(judge==3),1)>0
%     ansy2 = zeros(size(y(:,1)));
%     ansy2(find(judge==3)) = A(find(judge==3),:)*sys(:,2);
%     ansy2(ansy2==0)=nan;
%     plot(input,ansy2,'g.','MarkerSize',20)
% end
%
% if size(find(judge==4),1)>0
%     ansy3 = zeros(size(y(:,1)));
%     ansy3(find(judge==4)) = A(find(judge==4),:)*sys(:,3);
%     ansy3(ansy3==0)=nan;
%     plot(input,ansy3,'k.','MarkerSize',20)
% end
%
%
% legend('true y','Subysystem_1','Subysystem_2','Subysystem_3')
%


% %%
% tmp = x4;
% x4 = x3;
% x3 = tmp;
%%

Phi2 = [ones(size(y)) 1./x1 x2 sin(x1) cos(x2) exp(x1.*x2) x3  x4 x5];
para_log.idx_sys = idx_sys;
para_log.beta = 0.5;
para_log.y = y;
para_log.Phi2 = Phi2;
para_log.normalize = 1;
[syslogic,labelMat,data] = ihydelogic(para_log);

