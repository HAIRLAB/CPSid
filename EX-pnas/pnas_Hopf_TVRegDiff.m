% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

clear all, close all, clc
addpath('./utils');
addpath('./tools');
% generate Data
polyorder = 5;
usesine = 0;
eps = 0.005; 
lambda = 0.25;
mu = .1;
omega = 1;
A = 1;

dt = 0.0025;

tspan=[dt:dt:75];

N = length(tspan);
options = odeset('RelTol',1e-12,'AbsTol',1e-12*[1 1 ]);
x = [];
dx = [];
dxt = [];

%%  build up data over many values of bifurcation parameter
for mu = [-.15 -.05]  % start with stable mu
    clear tt xt dxt
    x0=[2; 0];  % Initial condition
    % Integrate
    [tt,xt]=ode45(@(t,x) hopf(t,x,mu,omega,A),tspan,x0,options);
    xtrue = xt;
    for i=1:length(xtrue)
        dxtrue(i,:) = hopf(0,xtrue(i,:),mu,omega,A);
    end
    xt = xt + eps*randn(size(xt));
    xt = [xt 0*xt(:,1)+mu];
    dxt(:,1) = TVRegDiff( xt(:,1), 5, 2, [], 'small', 1e+2, dt, 1, 1 );
    dxt(:,2) = TVRegDiff( xt(:,2), 5, 2, [], 'small', 1e+2, dt, 1, 1 );
    dxt = dxt(1:end-1,:);
    x = [x; xt(1000:end-500,:)];
    dxt(:,3) = 0*dxt(:,1);
    dx = [dx; dxt(1000:end-500,:)];
    plot(dxt(:,1:2),'k')
    hold on
    plot(dxtrue,'r')
    pause(1)
    hold off
end

for mu = [.05 .15 .25 .35 .45 .55]  % now use unstable mu
    % two ICs for this case... inside and outside limit cycle
    clear tt xt dxt
    x0=[.01; 0];  % Initial condition
    % Integrate
    [tt,xt]=ode45(@(t,x) hopf(t,x,mu,omega,A),tspan,x0,options);
    xtrue = xt;
    for i=1:length(xtrue)
        dxtrue(i,:) = hopf(0,xtrue(i,:),mu,omega,A);
    end    
    xt = xt + eps*randn(size(xt));
    xt = [xt 0*xt(:,1)+mu];
    dxt(:,1) = TVRegDiff( xt(:,1), 5, 10, [], 'small', 1e+1, dt, 1, 1 );
    dxt(:,2) = TVRegDiff( xt(:,2), 5, 10, [], 'small', 1e+1, dt, 1, 1 );
    dxt = dxt(1:end-1,:);
    x = [x; xt(1000:end-500,:)];
    plot(dxt(:,1:2),'k')
    hold on
    plot(dxtrue,'r')
    pause(1)
    hold off
    
    dxt(:,3) = 0*dxt(:,1);
    dx = [dx; dxt(1000:end-500,:)];
    clear tt xt dxt
    x0=[2; 0];  % Initial condition
    % Integrate
    [tt,xt]=ode45(@(t,x) hopf(t,x,mu,omega,A),tspan,x0,options);
    xtrue = xt;
    for i=1:length(xtrue)
        dxtrue(i,:) = hopf(0,xtrue(i,:),mu,omega,A);
    end    
    xt = xt + eps*randn(size(xt));
    xt = [xt 0*xt(:,1)+mu];
    dxt(:,1) = TVRegDiff( xt(:,1), 5, 10, [], 'small', 1e+1, dt, 1, 1 );
    dxt(:,2) = TVRegDiff( xt(:,2), 5, 10, [], 'small', 1e+1, dt, 1, 1 );
    dxt = dxt(1:end-1,:);
    x = [x; xt(1000:end-500,:)];
    plot(dxt(:,1:2),'k')
    hold on
    plot(dxtrue,'r')
    pause(1)
    hold off
    dxt(:,3) = 0*dxt(:,1);
    dx = [dx; dxt(1000:end-500,:)];
end
%%
% pool Data
Theta = poolData(x,3,polyorder,usesine);
m = size(Theta,2);

%%  iterative least squares solution
lambda = 0.85;
Xi0 = sparsifyDynamics(Theta,dx,lambda,3)
index = 1:100:size(dx,1);
Theta = Theta(index,:);
dx = dx(index,:);
parameter.lambda = [1 0.003];   % the lambda of z in algorithm 1.
parameter.MAXITER = 5;
parameter.max_s = 20;%the max s
parameter.epsilon = [1e-4  0.1];
parameter.Phi = Theta;
parameter.normalize_y = 1;
for i =1:3
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
poolDataLIST({'x','y','u'},Xi,3,polyorder,usesine);

%% Generate new data using sparse identified model
options = odeset('RelTol',1e-8,'AbsTol',1e-8*ones(1,3));

x = [];
for mu = [-.15 -.05]
    clear tt xt dxt     
    tspan=[dt:dt:75];
    x0=[2; 0; mu];  % Initial condition
    % Integrate
    [tt,xt]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),tspan,x0,options);
    x = [x; xt(:,:)];    
end

tspan=[dt:dt:75];
for mu = [.05 .15 .25 .35 .45 .55]
    clear tt xt dxt
    x0=[.01; 0; mu];  % Initial condition
    % Integrate
    [tt,xt]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),tspan,x0,options);
    x = [x; xt(:,:)];
    clear tt xt dxt
    x0=[2; 0; mu];  % Initial condition
    % Integrate
    [tt,xt]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),tspan,x0,options);
    x = [x; xt(:,:)];    
end


%% FIGURES
L = 30000;

Rv = 0:.01:sqrt(.6);
Tv = 0:2*pi/101:2*pi;
[R,T] = meshgrid(Rv,Tv);
X = R.*cos(T);
Y = R.*sin(T);
figure
zbottom = 0;

plot3([-.2 0],[0 0],[0 0],'k','LineWidth',2.5);hold on
plot3([0 .6],[0 0],[0 0],'k--','LineWidth',2.5);

h1=surf(X.^2+Y.^2,X,Y);
set(h1,'EdgeColor','none','FaceColor',[.5 .5 .5],'FaceAlpha',0.8)
lighting phong

tspan = 0:.01:70;
lambda = -100;
mu = 0.1;
omega = 1;
A = -mu;
k = 0.25;

% figure
for k=1:length(x)/L
    if(mod(k,2)&&k>2)
        plot3(x(1+(k-1)*L:k*L,3),x(1+(k-1)*L:k*L,1)*.999,x(1+(k-1)*L:k*L,2)*.999,'r','LineWidth',2)
    end
    if((~mod(k,2))||(k==1))
        plot3(x(1+(k-1)*L:k*L,3),x(1+(k-1)*L:k*L,1)*1.001,x(1+(k-1)*L:k*L,2)*1.001,'b','LineWidth',2)
    end
    hold on
    view(13,28)
end
view(10,28)

xlim([-.2 .6])
ylim([-1 1])
zlim([-1 1])
set(gca,'xtick',[-.2 0 .2 .4 .6],'xticklabel',{});
set(gca,'ytick',[-1 0 1],'yticklabel',{});
set(gca,'ztick',[-1 0 1],'zticklabel',{});
set(gca,'LineWidth',2)
grid on