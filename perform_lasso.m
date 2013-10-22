function a = perform_lasso(y,Phi,lambda,primal)

% perform_lasso - solve primal or dual lasso problems
%
%   a = perform_lasso(y,Phi,lambda,primal)
%
%   If primal=1, it solves the primal problem
%       min_a 1/2*|y-Phi*a|^2 + lambda*|a|_1  (lambda>0)
%       min_{y-Phi*a} |x|_1     (lambda=0)
%   If primal=-1, it solves the dual problem
%       min_a |y/lambda-a|^2  s.t. |Phi'*a|_inf<=1  (lambda>0)
%       min_a <y,a>   s.t.   |Phi'*a|_inf<=1  (lambda=0)
%
%   Copyright (c) 2013 Gabriel Peyre

if nargin<4
    primal = 1;
end

[p,n] = size(Phi);
    
if primal==1
    %% Primal problem %%    
    cvx_solver sdpt3 % SeDuMi %
    cvx_begin quiet
    cvx_precision high;
    variable a(n);
    if lambda==0
        Phi*a==y;
        minimize( norm(a,1) )
    else
%        minimize( 1/2*norm(Phi*a-y, 'fro')^2 + lambda*norm(a,1) )
        minimize( 1/2*sum(pow_abs(Phi*a-y,2)) + lambda*norm(a,1) )
    end
    cvx_end
else
    %% Dual Problem %%    
    cvx_solver sdpt3 % SeDuMi %
    cvx_begin quiet
    cvx_precision high;
    variable a(p) complex;
    norm(Phi'*a,Inf)<=1;
    if lambda==0
        maximize(real( a' * y ))
    else
        minimize( norm(y/lambda-a, 'fro') )
    end
    cvx_end
end
