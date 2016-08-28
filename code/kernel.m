function y = kernel(name, t, k, options)

% kernel - Load a low-frequency kernel
%
%   y = kernel(name, t, k, options);
%
%   If name='dirichlet', compute
%       y = d^k K_{n,a}(t)/dt^k
%   where
%       K_{n,lambda}(t) = sin(a*t)^n / sin(pi*t)^n
%
%   Compute, for t in [0,1]
%       sum_{f=-fc}^{fc} (2 i pi f)^(k) w(f+fc+1) exp( 2*i*pi * fc * t )
%   so that k is the order of the derivative.
%
%   Copyright (c) 2013 Gabriel Peyre

options.null = 0;
fc = getoptions(options, 'fc', -1);
if fc<0
    warning('You need to set options.fc');
end

switch name
        
    case 'dirichlet'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        n = getoptions(options, 'power', 1);
        a = getoptions(options, 'a', -1);
        if a<0
            a = (2*fc/n+1)*pi;
        end

        if k==2
            t = t + 1.12345678*1e-6; % avoid division by zero
        end
        
        ca = cos(a*t); sa = sin(a*t);
        cp = cos(pi*t); sp = sin(pi*t);
        
        switch k
            case 0
                y = sin(a*t).^n ./ sin(pi*t).^n;
                y(t==0)=(a/pi)^n;
            case 1
                options.a = a;
                options.power = n-1;
                u = kernel(name, t, 0, options) ./ sp.^2;
                v = a*ca.*sp - pi*sa.*cp;
                %
                y = n * u.*v;
                y(t==0)=0;
            case 2
                options.a = a;
                options.power = n-1;
                u = kernel(name, t, 0, options) ./ sp.^2;
                u1 = ( kernel(name, t, 1, options).*sp.^2 - ...
                    2*pi*cp.*sp.*kernel(name, t, 0, options) ) ./ sp.^4;
                %
                v = a*ca.*sp - pi*sa.*cp;
                v1 = - (a^2-pi^2) * sa .* sp;
                %
                y = n * (u1.*v + u.*v1);
            otherwise
                error('Only derivatives up to order 2 are implemented.');
        end
        
    
    
    case 'gaussian'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sigma = 1/(4*fc);
        t(t>1/2) = t(t>1/2)-1;
        switch k
            case 0
                y = exp( -t.^2/(2*sigma^2) );
            case 1
                y = -t/sigma^2 .* exp( -t.^2/(2*sigma^2) );
            case 2
                y = (t.^2-1)/sigma^2 .* exp( -t.^2/(2*sigma^2) );
            otherwise
                error('Only derivatives up to order 2 are implemented.');
        end
        
        
    case 'fourier'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        w = getoptions(options, 'w', ones(2*fc+1,1));
        w = w(:)';        
        N = length(t(:));
        omega = -fc:fc;
        w1 = w .* (2i*pi.*omega).^k;
        [Omega,U] = meshgrid(omega,t(:));
        Z = exp(2i*pi*Omega.*U) .* repmat(w1, [N 1]);
        y = reshape( real(sum(Z,2)), size(t) );
        
    otherwise
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        error('Unknown');
        
end
