function w_estimate = sparsesolver(Output, Dic, lambda, MAXITER, Omega,ifprint)
%%
% This is a VANILA implementation of the follwing paper
% 
% The algorithm solves the inverse problem 
%   $y = Aw + \xi$
%
% ============= Inputs ===============
% y                                    : output, in the paper $y(t) =
%                                           (x(t+\delta)-x(t))/\delta$;
% A                                    : dictionary matrix;
% lambda                         : the tradeoff parameter you should use, 
%                                          basically, it is proportional to
%                                          the invese of the variance, e.g. 1;
% MAXITER                    : maximum number of the iterations, e.g. 5

% ============= Reference =============
% W. Pan, Y. Yuan, J. Goncalves, and G.-B. Stan, 
% A Sparse Bayesian Approach to the Iden- tification of Nonlinear State-Space Systems,
% IEEE Transaction on Automatic Control, 2015 (to appear). arXiv:1408.3549
% http://arxiv.org/abs/1408.3549
%
% ============= Author =============
%  Wei Pan (w.pan11@imperial.ac.uk, panweihit@gmail.com)
%
% ============= Version =============
%   1.0 (Sep ?, 2012)
%%

delta = 1e-5;

[M,N]=size(Dic);
% initialisation of the variables
U=ones(N, MAXITER);
Gamma=zeros(N, MAXITER);
UU=zeros(N, MAXITER);
w_estimate=zeros(N, MAXITER);
WWW=ones(N, MAXITER);
if ifprint ==1
    fprintf(1, 'Finding a sparse feasible point using l1-norm heuristic ...');
end

for iter=1:1:MAXITER
    if ifprint==1
        fprintf('This is round %d \n', iter);
    end
    cvx_begin quiet
    cvx_solver sedumi   %sdpt3
    variable W(N)
    minimize    (lambda*norm( U(:,iter).*W, 1)+ 0.5*(Dic*W-Output)'*Omega*(Dic*W-Output))
    %                 subject to
    %                           W.^2-ones(101,1)<=0;
    cvx_end
    
    if isnan(W)
        W = w_estimate(:,iter-1);
        
    end
    w_estimate(:,iter)=W;
    WWW(:,iter)=W;
    Gamma(:,iter)=U(:,iter).^-1.*abs(W);
    Dic0=lambda*eye(M)+Dic*diag(Gamma(:,iter))*Dic';
    UU(:, iter)=diag(Dic'*(Dic0\Dic));
    U(:,iter+1)=abs(sqrt(UU(:, iter)));
    if isnan(U(:,iter+1))
        U(:,iter+1) = U(:,iter);
        
    end

    
    
%     for i=1:N
%         if   w_estimate(i,iter).^2/norm(w_estimate(:,iter))^2<delta
%             w_estimate(i,iter)=0;
%         end
%     end
     
end

