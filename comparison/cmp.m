function theta = ref18(parameter)

y = parameter.y;
Phi = parameter.Phi;
lambda = parameter.lambda;
[M,N] = size(Phi);

cvx_begin quiet
cvx_solver sedumi
variable theta(N)
variable e(M)
minimize (lambda*norm((y-Phi*theta-e),1)+0.5*norm(e,2))
cvx_end

end