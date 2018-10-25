function dy = hopf(t,y,mu,omega,A)
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

dy = [
mu*y(1) - omega*y(2) - A*y(1)*(y(1)^2+y(2)^2);
omega*y(1) + mu*y(2) - A*y(2)*(y(1)^2+y(2)^2);
];