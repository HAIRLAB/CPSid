% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

clear all, close all, clc
addpath('./utils');
figpath = '../figures/';
addpath('../tools');
%% generate Data
polyorder = 5;
usesine = 0;
n = 3;  % 3D system
A = [-.1 2 0; -2 -.1 0 ; 0 0 -.3];
rhs = @(x)A*x;   % ODE right hand side
tspan=[0:.01:10];   % time span
x0 = [2; 0; 1];        % initial conditions
options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n));
[t,x]=ode45(@(t,x)rhs(x),tspan,x0,options);  % integrate

%% compute Derivative
eps = .01;
for i=1:length(x)
    dx(i,:) = A*x(i,:)';
end
dx = dx + eps*randn(size(dx));

%% pool Data  (i.e., build library of nonlinear time series)
Theta = poolData(x,n,polyorder,usesine);
m = size(Theta,2);

%% compute Sparse regression: sequential least squares
parameter.lambda = [1 0.002];   % the lambda of z in algorithm 1.
parameter.MAXITER = 5;
parameter.max_s = 20;%the max s
parameter.epsilon = [1e-4  0.05];

parameter.Phi = Theta;

parameter.normalize_y = 1;


for i =1:size(dx,2)
    parameter.y = dx(:,i);
    [result]=ihyde(parameter);
    
    if size(result.sys,2)>1
        result.epsilon = parameter.epsilon(2);
        result.lambda = parameter.lambda(2);
        result.threshold = [0.05];
        final_result  = finetuning( result);

    else 
        final_result = result;
    end
    sys = final_result.sys;
    idx_sys = final_result.idx;
    Xi(:,i) = sys;
    
end

poolDataLIST({'x','y','z'},Xi,n,polyorder,usesine);

%% integrate true and identified systems
[tA,xA]=ode45(@(t,x)rhs(x),tspan,x0,options);   % true model
[tB,xB]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),tspan,x0,options);  % approximate

%% FIGURES!!
figure
dtA = [0; diff(tA)];
plot3(xA(:,1),xA(:,2),xA(:,3),'r','LineWidth',1.5);
view(49,20)
hold on
dtB = [0; diff(tB)];
plot3(xB(1:5:end,1),xB(1:5:end,2),xB(1:5:end,3),'k--','LineWidth',1.2);
l1 = legend('True','Identified');
set(l1,'FontSize',13,'Location','NorthEast');
set(gca,'FontSize',13);
xlabel('x_1','FontSize',13)
ylabel('x_2','FontSize',13)
zlabel('x_3','FontSize',13)
view(49,20)
grid on

figure
plot(tA,xA(:,1),'r','LineWidth',1.5)
hold on
plot(tA,xA(:,2),'b','LineWidth',1.5)
plot(tA,xA(:,3),'g','LineWidth',1.5)
plot(tB,xB(:,1),'k--','LineWidth',1.2)
plot(tB,xB(:,2),'k--','LineWidth',1.2)
plot(tB,xB(:,3),'k--','LineWidth',1.2)
xlabel('Time','FontSize',13)
ylabel('State, x_k','FontSize',13)
legend('True x_1','True x_2','True x_3','Identified')
