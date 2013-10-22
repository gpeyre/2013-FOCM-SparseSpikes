function p = perform_sdp_dual(y,lambda, options)

% perform_sdp_dual - compute dual vector
%
%   p = perform_sdp_dual(y,lambda,options);
%
%   If lambda=0, solves
%       min_{p} <p,y>   s.t.  |F^* diag(conj(h)) p|_inf<=1
%   If lambda>0, solves
%       min_{p} <p,y>   s.t.  |F^* diag(conj(h)) p|_inf<=1 
%   where F^* is the inverse Fourier transform 
%       F^* p(t) = sum_{f=-fc}^fc p(f) e^{2*i*pi*f*t}
%
%   options.h : fft of the filter
%   options.solver : either 'cvx' or 'dr' (Douglas-Rachford splitting).
%
%   The idea of solving the dual problem using a SDP lifting is introduced
%   in: 
%       Towards a Mathematical Theory of Super-Resolution
%       Emmanuel J. Candes and Carlos Fernandez-Granda
%
%   Copyright (c) 2013 Gabriel Peyre

n = length(y);

options.null = 0;
h = getoptions(options, 'h', ones(n,1));
solver = getoptions(options, 'solver', 'cvx');

switch solver
    case 'cvx'
        cvx_solver sdpt3 % SeDuMi %
        cvx_begin sdp quiet
        cvx_precision high;
        variable X(n+1,n+1) hermitian;
        variable p(n) complex;
        X >= 0;
        X(n+1,n+1) == 1;
        X(1:n,n+1) == p .* conj(h);
        trace(X) == 2;
        for j = 1:n-1,
            sum(diag(X,j)) == X(n+1-j,n+1);
        end
        if lambda==0
            maximize(real( p' * y ))
        else
            minimize( norm(y/lambda-p, 'fro') )
        end
        cvx_end
        
    case 'dr'
        %%% DOUGLAS RACHFORD %%%
        % min F(X) + G(X)
        %   where
        %       F(X) = f(p(X)) + i_{SOS}(Q(X))
        %       G(X) = i_{SDP}(X)
        %   where
        %       X = [Q, p(Q); p(Q)', 1];
        % helper
        dotp = @(x,y)real(x'*y);
        Xmat = @(p,Q)[Q, p; p', 1];
        Qmat = @(X)X(1:end-1,1:end-1);
        pVec = @(X)X(1:end-1,end);
        % function to minimize
        if lambda>0
            f = @(p)1/2*norm( y/lambda-p )^2;
            Proxf = @(p,gamma)( p + gamma*y/lambda )/(1+gamma);
        else
            f = @(p)-dotp(y,p);
            Proxf = @(p,gamma)p + gamma*y;
        end
        % prox operators
        ProxF = @(X,gamma)Xmat( Proxf(pVec(X),gamma/2), perform_sos_projection(Qmat(X)) );
        ProxG = @(X,gamma)perform_sdp_projection(X);
        % initial point
        X0 = zeros(n+1);
        options.niter = getoptions(options, 'niter', 500);
        options.gamma = getoptions(options, 'gamma', 1/10);
        options.verb = getoptions(options, 'verb', 0);
        debug_dr = getoptions(options, 'debug_dr', 0);
        if debug_dr
            options.report = @(X)struct('ObjVal', f(pVec(X)), ...
                'ConstrSDP', min(real(eig(X))), ...
                'ConstrSOS', norm(perform_sos_projection(Qmat(X))-Qmat(X), 'fro') ); 
                % 'Error', norm(pVec(X)-p0),  ...
        end
        [X,R] = perform_dr(X0,ProxF,ProxG,options);
        p = pVec(X);
        if debug_dr
            ObjVal = s2v(R, 'ObjVal');
            ConstrSDP = s2v(R, 'ConstrSDP');
            ConstrSOS = s2v(R, 'ConstrSOS');
            % Error = s2v(R, 'Error');
            sel = 4:options.niter;
            clf;
            subplot(3,1,1);
            plot(ObjVal(sel)); axis tight; title('Objective');
            subplot(3,1,2);
            plot(ConstrSDP(sel)); axis tight; title('SDP');
            subplot(3,1,3);
            plot(ConstrSOS(sel)); axis tight; title('SOS');
        end
        
    otherwise
        error('Unknown solver');
        
end
        
end

function X = perform_sdp_projection(X)

% perform_sdp_projection - project on the cone of SDP matrices.
%
%   X = perform_sdp_projection(X);
%
%   Copyright (c) 2013 Gabriel Peyre

[V,D] = eig(X);
X = V*max(real(D),0)*V';

end

function Q = perform_sos_projection(Q)

% perform_sos_projection - project on the linear constraint for SOS polynomials
%
%   Q = perform_sos_projection(Q);
%
%   The constraint is that each non-leading diagonal should sum to 0 and
%   the main diagonal should sum to 0.
%
%   Copyright (c) 2013 Gabriel Peyre

n = size(Q,1);
% index set
U = reshape(1:n^2, [n n]);
for i=-n+1:n+1
    I = diag(U,i);
    Q(I) = Q(I) - mean(Q(I));
end
% add back id/n so that diag sum to 1.
Q = Q + eye(n)/n;

end
