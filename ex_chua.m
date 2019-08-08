clear all,close all,clc

addpath('./tools');
addpath('./data');
addpath('./SLR_dev');
%%  generate Data

data=load('chua.mat');
dt = (data.s(2,1) - data.s(1,1));

win_size = 30;
data.CH1V = wiener2(data.CH1V,[win_size 1]);
data.CH2V = wiener2(data.CH2V,[win_size 1]);
N=50;

y1 = data.CH1V(mod(0:length(data.CH1V)-1,1+N)<1);
y2 = data.CH2V(mod(0:length(data.CH2V)-1,1+N)<1);

t = 1:size(y1,1)*1.02e-7;
[y1, dy1, index] = estimatediff(y1, t, 'solver', 1, []);
y2 = y2(index);

index = 31:207;%because of wiener
y1 = y1(index);
y2 = y2(index);

dy1 = dy1(index);

ans_sys_idx{1} = [];
ans_sys_idx{2} = [];
ans_sys_idx{3} = [];
for i = 1:size(y1,1)
    if y1(i,1)<-1.5&&y1(i,1)>-6
        ans_sys_idx{1} = [ans_sys_idx{1} i ];
        state(i) = 1;
    end
    if y1(i,1)<1.5&&y1(i,1)>-1.5
        ans_sys_idx{2} = [ans_sys_idx{2} i ];
        state(i) = 2;
    end
    if y1(i,1)<6&&y1(i,1)>1.5
        ans_sys_idx{3} = [ans_sys_idx{3} i ];
        state(i) = 3;
    end
end


%%


A = [ones(size(y1)) y1 y2 exp(y1) y1./y2  cos(0.1*y1).^2./(1+y2.^2)  cos(y1+y2).^2];


parameter.lambda = [0.05 0.01];    % the lambda of z in algorithm 1.
parameter.MAXITER = 5;
parameter.max_s = 20;%the max s
parameter.epsilon = [0.012  0.044];
parameter.Phi = A;
parameter.y = dy1;
parameter.normalize_y = 1;
[result]=ihyde(parameter);



result.epsilon = parameter.epsilon(2);
result.lambda = parameter.lambda(2);
result.threshold = [0.05];
final_result  = finetuning( result);
sys = final_result.sys;
idx_sys = final_result.idx;

%%
Phi2 = [ ones(size(y1)) y1 sin(y2).*cos(y1)  dy1./(sin(y1)+dy1) dy1./y2  dy1];

Phi2 = [ ones(size(y1)) y1 ];

para_log.idx_sys = idx_sys;
para_log.beta = 0.01;

para_log.y = dy1;
para_log.Phi2 = Phi2;

para_log.normalize = 1;

[syslogic,labelMat,data] = ihydelogic(para_log);

for i =1:size(para_log.idx_sys,2)
    for j=1:size(para_log.idx_sys,2)
        if length(syslogic{i,j})~=0
            logicsys(i,j) = -syslogic{i,j}(1)/syslogic{i,j}(2);
        end
    end
    
end

%%
judge = zeros(size(state));
judge(idx_sys{1}) =3;
judge(idx_sys{2}) =1;
judge(idx_sys{3}) =2;
wrong_numbers = find((judge - state)~=0);
%%
 close all
close all
figure;

% Create axes
axes1 = axes;
hold(axes1,'on');
index1 = 1:177;
y11 = dy1(index1);
judge11 = state(index1);



ansy1 = 1000000*ones(size(y11(:,1)));
ansy1(find(judge11==1)) = y11(find(judge11==1));
ansy1(ansy1==1000000)=nan;
plot(ansy1,'MarkerSize',4,'LineWidth',3,...
    'Color',[0.925490200519562 0.839215695858002 0.839215695858002]);

ansy11 = 1000000*ones(size(y11(:,1)));
ansy11(idx_sys{2}) = A(idx_sys{2},:)*sys(:,2);
ansy11(ansy11==1000000)=nan;
plot(ansy11,'MarkerSize',8,'Marker','o','LineWidth',1.5,'LineStyle','none',...
    'Color',[0.674509823322296 0.164705887436867 0.0392156876623631])
err1 = ansy1 - ansy11;

ansy2 = 1000000*ones(size(y11(:,1)));
ansy2(find(judge11==2)) = y11(find(judge11==2));
ansy2(ansy2==1000000)=nan;
plot(ansy2,'MarkerSize',4,'LineWidth',3,...
    'Color',[0.600000023841858 0.800000011920929 1])

ansy22 = 1000000*ones(size(y11(:,1)));
ansy22(idx_sys{3}) = A(idx_sys{3},:)*sys(:,3);
ansy22(ansy22==1000000)=nan;
plot(ansy22,'MarkerSize',8,'Marker','o','LineWidth',1.5,'LineStyle','none',...
    'Color',[0.223529413342476 0.396078437566757 0.737254917621613])
err2 = ansy2 - ansy22;


ansy3 = 1000000*ones(size(y11(:,1)));
ansy3(find(judge11==3)) = y11(find(judge11==3));
ansy3(ansy3==1000000)=nan;
plot(ansy3,'MarkerSize',4,'LineWidth',3,...
    'Color',[0.756862759590149 0.866666674613953 0.776470601558685])

ansy33 = 1000000*ones(size(y11(:,1)));
ansy33(idx_sys{1}) = A(idx_sys{1},:)*sys(:,1);
ansy33(ansy33==1000000)=nan;
plot(ansy33,'MarkerSize',8,'Marker','o','LineWidth',1.5,'LineStyle','none',...
    'Color',[0.301960796117783 0.686274528503418 0.290196090936661])
err3 = ansy3 - ansy33;
% Create xlabel
xlabel({'t'});

% Create ylabel
ylabel({'$dy_1$'},'Interpreter','latex');

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',24,'LineWidth',2.5);
