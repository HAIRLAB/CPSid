% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

% Note, for larger error terms, it helps to remove constant terms in
% "poolData"

clear all, close all, clc
figpath = '../figures/';
addpath('./utils');
addpath('./tools');
%% generate Data
polyorder = 5;
usesine = 0;
sigma = 10;  % Lorenz's parameters (chaotic)
beta = 8/3;
rho = 28;
n = 3;
x0=[-8; 8; 27];  % Initial condition

% Integrate
dt = 0.001;
tspan=[dt:dt:50];
N = length(tspan);
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
[t,x]=ode45(@(t,x) lorenz(t,x,sigma,beta,rho),tspan,x0,options);
xclean = x;
% add noise
eps = .01;
x = x + eps*randn(size(x));
% compute clean derivative  (just for comparison!)
for i=1:length(x)
    dxclean(i,:) = lorenz(0,xclean(i,:),sigma,beta,rho);
end

%%  Total Variation Regularized Differentiation
dxt(:,1) = TVRegDiff( x(:,1), 10, .00002, [], 'small', 1e12, dt, 1, 1 );
hold on
plot(dxclean(:,1),'r')
xlim([5000 7500])
figure
dxt(:,2) = TVRegDiff( x(:,2), 10, .00002, [], 'small', 1e12, dt, 1, 1 );
hold on
plot(dxclean(:,2),'r')
xlim([5000 7500])
figure
dxt(:,3) = TVRegDiff( x(:,3), 10, .00002, [], 'small', 1e12, dt, 1, 1 );
hold on
plot(dxclean(:,3),'r')
xlim([5000 7500])

%
xt(:,1) = cumsum(dxt(:,1))*dt;
xt(:,2) = cumsum(dxt(:,2))*dt;
xt(:,3) = cumsum(dxt(:,3))*dt;
xt(:,1) = xt(:,1) - (mean(xt(1000:end-1000,1)) - mean(x(1000:end-1000,1)));
xt(:,2) = xt(:,2) - (mean(xt(1000:end-1000,2)) - mean(x(1000:end-1000,2)));
xt(:,3) = xt(:,3) - (mean(xt(1000:end-1000,3)) - mean(x(1000:end-1000,3)));
xt = xt(1000:end-1000,:);
dxt = dxt(1000:end-1000,:);  % trim off ends (overly conservative)

%% pool Data  (i.e., build library of nonlinear time series)
Theta = poolData(xt,n,polyorder,usesine);
m = size(Theta,2);

%% compute Sparse regression: sequential least squares
lambda = 0.02;      % lambda is our sparsification knob.
index = 1:25:length(dxt);
Theta = Theta(index,:);
dxt = dxt(index,:);

Xi0 = sparsifyDynamics(Theta,dxt,lambda,n)
%%
parameter.lambda = [1 0.00365];   % the lambda of z in algorithm 1.
parameter.MAXITER = 5;
parameter.max_s = 20;%the max s
parameter.epsilon = [1e-4  2];
parameter.Phi = Theta;
parameter.normalize_y = 1;
for i =1:size(dxt,2)
    parameter.y = dxt(:,i);
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
%% FIGURE 1:  LORENZ for T\in[0,20]
tspan = [0 20];
[tA,xA]=ode45(@(t,x)lorenz(t,x,sigma,beta,rho),tspan,x0,options);   % true model
[tB,xB]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),tspan,x0,options);  % approximate

figure
subplot(1,2,1)
dtA = [0; diff(tA)];
color_line3(xA(:,1),xA(:,2),xA(:,3),dtA,'LineWidth',1.5);
view(27,16)
grid on
xlabel('x','FontSize',13)
ylabel('y','FontSize',13)
zlabel('z','FontSize',13)
set(gca,'FontSize',13)
subplot(1,2,2)
dtB = [0; diff(tB)];
color_line3(xB(:,1),xB(:,2),xB(:,3),dtB,'LineWidth',1.5);
view(27,16)
grid on

% Lorenz for t=20, dynamo view
figure
subplot(1,2,1)
plot(tA,xA(:,1),'k','LineWidth',1.5), hold on
plot(tB,xB(:,1),'r--','LineWidth',1.5)
grid on
xlabel('Time','FontSize',13)
ylabel('x','FontSize',13)
set(gca,'FontSize',13)
subplot(1,2,2)
plot(tA,xA(:,2),'k','LineWidth',1.5), hold on
plot(tB,xB(:,2),'r--','LineWidth',1.5)
grid on


%% FIGURE 1:  LORENZ for T\in[0,250]
tspan = [0 250];
options = odeset('RelTol',1e-6,'AbsTol',1e-6*ones(1,n));
[tA,xA]=ode45(@(t,x)lorenz(t,x,sigma,beta,rho),tspan,x0,options);   % true model
[tB,xB]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),tspan,x0,options);  % approximate

figure
subplot(1,2,1)
dtA = [0; diff(tA)];
color_line3(xA(:,1),xA(:,2),xA(:,3),dtA,'LineWidth',1.5);
view(27,16)
grid on
xlabel('x','FontSize',13)
ylabel('y','FontSize',13)
zlabel('z','FontSize',13)

subplot(1,2,2)
dtB = [0; diff(tB)];
color_line3(xB(:,1),xB(:,2),xB(:,3),dtB,'LineWidth',1.5);
view(27,16)
grid on
xlabel('x')
ylabel('y')
zlabel('z')