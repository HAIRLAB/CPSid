% This file works with MATLAB and is automatically generated with
% the System Biology Format Converter (http://sbfc.sourceforge.net/)
% from an SBML file.
% To run this file with Octave you must edit the comments providing
% the definition of the ode solver and the signature for the
% xdot function.
%
% The conversion system has the following limitations:
%  - You may have to re order some reactions and Assignment Rules definition
%  - Delays are not taken into account
%  - You should change the lsode parameters (start, end, steps) to get better results
%

%
% Model name = Courtemanche1998_AtrialActionPotential
%
% is http://identifiers.org/biomodels.db/MODEL0913049417
% isDescribedBy http://identifiers.org/pubmed/9688927
%
clear all
close all
addpath('../tools');
addpath('../data');


%Initial conditions vector
x0=zeros(21,1);
x0(1) = -81.18;
x0(2) = 0.002908;
x0(3) = 0.9649;
x0(4) = 0.9775;
x0(5) = 0.03043;
x0(6) = 0.9992;
x0(7) = 0.004966;
x0(8) = 0.9986;
x0(9) = 3.296E-5;
x0(10) = 0.01869;
x0(11) = 1.367E-4;
x0(12) = 0.9996;
x0(13) = 0.7755;
x0(14) = 2.35E-112;
x0(15) = 1.0;
x0(16) = 0.9992;
x0(17) = 11.17;
x0(18) = 139.0;
x0(19) = 1.013E-4;
x0(20) = 1.488;
x0(21) = 1.488;


% Depending on whether you are using Octave or Matlab,
% you should comment / uncomment one of the following blocks.
% This should also be done for the definition of the function f below.
% Start Matlab code
tspan=[0:0.005:600];
% opts = odeset('AbsTol',1e-3,'Maxstep',0.01);
[t,x]=ode45(@AP_simulation,tspan,x0);
% plot(t,x);
% End Matlab code

% Start Octave code
%	t=linspace(0,100,100);
%	x=lsode('f',x0,t);
% End Octave code
for i = 1:size(x,2)
    [x_tmp1, x_tmp2, index] = estimatediff(x(:,i), t, 'solver', 1, []);
    x_dot(:,i) = x_tmp2;
    x_tmp(:,i) = x_tmp1;
end
x = x(index,:);
t = t(index,:);
index = 24000:60:100000;
x_dot = x_dot(index,:);
x_tmp = x_tmp(index,:);
x = x(index,:);
t = t(index,:);

% plot(t,x);
x = x_tmp;
tmp = x(:,1);
for i = 1:size(x,1)
    global_par_alpha_h=piecewise(0.135*exp((tmp(i)+80)/(-6.8)), (40+tmp(i)) < 0, 0);
    global_par_beta_h=piecewise(3.56*exp(0.079*tmp(i))+310000*exp(0.35*tmp(i)), (40+tmp(i)) < 0, 1/(0.13*(1+exp((tmp(i)+10.66)/(-11.1)))));
    global_par_tau_h=1/(global_par_alpha_h+global_par_beta_h);
    global_par_h_inf=global_par_alpha_h/(global_par_alpha_h+global_par_beta_h);
    y3(i,1) =     (global_par_h_inf-x(i,3))/global_par_tau_h;
    global_par_alpha_j=piecewise(((-127140)*exp(0.2444*tmp(i))-3.474E-5*exp((-0.04391)*tmp(i)))*(tmp(i)+37.78)/(1+exp(0.311*(tmp(i)+79.23))), 40+tmp(i) < 0, 0);
    global_par_beta_j=piecewise(0.1212*exp((-0.01052)*tmp(i))/(1+exp((-0.1378)*(tmp(i)+40.14))), 40+tmp(i) < 0, 0.3*exp((-2.535E-7)*tmp(i))/(1+exp((-0.1)*(tmp(i)+32))));
    global_par_j_inf=global_par_alpha_j/(global_par_alpha_j+global_par_beta_j);
    global_par_tau_j(i,1)=1/(global_par_alpha_j+global_par_beta_j);
    y4(i,1) = (global_par_j_inf-x(i,4))/global_par_tau_j(i,1);
end

global_par_alpha_h1 = 0.135*exp((tmp+80)/(-6.8));
global_par_alpha_h2 = zeros(size(y3));
global_par_beta_h1 = 3.56*exp(0.079*tmp)+310000*exp(0.35*tmp);
global_par_beta_h2 = 1./(0.13*(1+exp((tmp+10.66)/(-11.1))));
global_par_alpha_j1=((-127140)*exp(0.2444*tmp)-3.474E-5*exp((-0.04391)*tmp)).*(tmp+37.78)./(1+exp(0.311*(tmp+79.23)));
global_par_alpha_j2 = zeros(size(y3));
global_par_beta_j1 = 0.1212*exp((-0.01052)*tmp)./(1+exp((-0.1378)*(tmp+40.14)));
global_par_beta_j2 = 0.3*exp((-2.535E-7)*tmp)./(1+exp((-0.1)*(tmp+32)));

% dictionary matrix of gating variable h 
A3 = [exp((tmp+80)/(-6.8)) exp((tmp+80)/(-6.8)).*x(:,3) (global_par_beta_h1).*x(:,3)  global_par_beta_h2.*x(:,3)];
% dictionary matrix of gating variable j 
% A4 = [global_par_alpha_j1 global_par_alpha_j1.*x(:,4)  (exp((-0.01052)*tmp)./(1+exp((-0.1378)*(tmp+40.14)))).*x(:,4)  (exp((-2.535E-7)*tmp)./(1+exp((-0.1)*(tmp+32)))).*x(:,4)];
polyorder = 3;
memory = 0;
usesine  = 0;
basis_function.work = 'off';
A= library([tmp,x(:,4)],polyorder,usesine,memory,basis_function);
state = zeros(size(y4));
state((40+tmp) < 0) = 1;
state((40+tmp) > 0) = 2;
parameter.lambda = [1e-4 3e-5];   % x3
parameter.lambda = [1e-4 3e-8];   % x3p5
% parameter.lambda = [1e-4 1e-5];   % x4
parameter.lambda = [1e-5 1e-4];   % x4p3
parameter.MAXITER = 5;
parameter.max_s = 20;   %the max s
parameter.epsilon = [3e-5  5e-5]; % for x(4) and x(3)
parameter.state = state;
parameter.Phi = A;
parameter.y = x_dot(:,4);
parameter.normalize_y = 1;
[result]=ihyde(parameter);
sys = result.sys;

%% finetuning 
result.epsilon = parameter.epsilon(2);
result.lambda = parameter.lambda(2);
result.threshold = [0.05];
final_result  = finetuning( result);
sys = final_result.sys;
idx_sys = final_result.idx;


% Depending on whether you are using Octave or Matlab,
% you should comment / uncomment one of the following blocks.
% This should also be done for the definition of the function f below.
% Start Matlab code

% adding few functions representing operators used in SBML but not present directly
% in either matlab or octave.
%% identify the transit logic 
Phi2 = [ones(size(x(:,1))) ,x(:,1)];

para_log.idx_sys = idx_sys;
para_log.beta = 1e-6;

para_log.y = y4;
para_log.Phi2 = Phi2;

para_log.normalize = 1;

[syslogic,labelMat,data] = ihydelogic(para_log);

%%  plot figure
subsys2 = find(state-2);
changetime = subsys2(1);
subsys2 = find(state-2);
changetime = subsys2(1);
figure()
subplot(2,1,1)
plot(t(2:changetime),x(1:changetime-1,3),'r','LineWidth',2);
hold on
plot(t(changetime+1:end),x(changetime:end-1,3),'b','LineWidth',2)
set(gca,'FontSize',12)
legend('subsystem_1','subsystem_2');
ylabel('h');
%text(49.6,0.38,'\leftarrow V>=-40','FontSize',14)
subplot(2,1,2)
plot(t(2:changetime),x(1:changetime-1,4),'r','LineWidth',2)
hold on 
plot(t(changetime+1:end),x(changetime:end-1,4),'b','LineWidth',2)
set(gca,'FontSize',12)
%text(49.60,0.64,'\leftarrow V>=-40','FontSize',14)
legend('subsystem_1','subsystem_2');
ylabel('j')
xlabel('time(ms)')
