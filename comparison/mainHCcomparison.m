clear
clc
close all

addpath('../tools');
addpath('../data');


%% design
fileID = fopen('hc_n6.txt');
formatSpec = '%s';
N = 3;
num = 500;
C_text = textscan(fileID,formatSpec,N,'Delimiter','|');
C_data0 = textscan(fileID,'%f %f %f ');
y = C_data0{1}(1:num,:);
x = C_data0{2}(1:num,:);
state = C_data0{3}(1:num,:);


%% make dic
A = [ones(size(x)) x x.^2  x.^3 x.^4 x.^5 ];


%% identify the subsystems 
% Bako's method
parameter1.lambda = 0.015; 
parameter1.epsilon = 0.1; 
parameter1.MAXITER = 1;
parameter1.max_s = 20;%the max s
parameter1.state = state;
parameter1.Phi = A;
parameter1.y = y;
parameter1.normalize_y = 1; 
result1 = Method1(parameter1) ;

% IHYDE with 1 iteration 
parameter2 = parameter1;
parameter2.lambda = [0.015 0.015];  
parameter2.epsilon = [0.25 0.2]; 
parameter2.MAXITER = 1;
result2 = ihyde(parameter2) ;

% IHYDE with 5 iteration 
parameter3 = parameter2;
parameter3.lambda = [0.03 0.008]; 
parameter3.epsilon = [1e-4 0.1]; 
parameter3.MAXITER = 5;
result3 = ihyde(parameter3) ;

%% comparsion results 
% display the identified systems
fprintf('The identified systems of Bako Method: \n')
display(cell2mat(result1.sys) )  
fprintf('The identified systems of IHYDE (1 iteration): \n')
display((result2.sys))
fprintf('The identified systems of IHYDE (5 iterations): \n')
display((result3.sys))

 
idx2 = find(state==2) ; 
idx3 = find(state==3) ; 
CorrectPoints1 = length(intersect(result1.idx{1},idx3))+length(intersect(result1.idx{2},idx2)) ;
CorrectPoints2 = length(intersect(result2.idx_sys{1},idx3))+length(intersect(result2.idx_sys{2},idx2)) ;
CorrectPoints3 = length(intersect(result3.idx_sys{1},idx3))+length(intersect(result3.idx_sys{2},idx2)) ;
fprintf('The number of misclassified points for Bako Method is %d \n',num-CorrectPoints1)
fprintf('The number of misclassified points for ihyde (1 iteration) is %d \n',num-CorrectPoints2)
fprintf('The number of misclassified points for ihyde (5 iterations) is %d \n',num-CorrectPoints3)

