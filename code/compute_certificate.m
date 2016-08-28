function [eta,etaF] = compute_certificate(w, x, s, options)

% compute_certificate - compute dual (pre-)certificate
%
%   [eta,etaF] = compute_certificate(phi, x, s, options);
%
%   eta is the Dual certificate interpolated on p points
%   etaF are the 2*fc+1 Fourer coefficients of eta (so this gives access to
%       the exact certificate)
%
%   w is the Fourier transform of the filter, of size 2*fc+1 
%   x stores the position of the Dirac's
%   s srotes the signs of the Dirac's
%
%   options.P = number of sampling points
%   options.order = number of vanishing derivatives
%   options.type = 'leastsquare' (for the least square certificate)
%                = 'fejer' (for the one developped by Candes and
%                Fernandez-Granda)
%                = 'positive' (only for positive s)
%
%   Copyright (c) 2013

options.null = 0;
order = getoptions(options, 'order', 1);
Type = getoptions(options, 'type', 'leastsquare');
P = getoptions(options, 'P', 2048);

x = x(:);
fc = (length(w)-1)/2;
n = length(x);
% observation grid
u = (0:P-1)'/P;

switch Type
    case 'leastsquare'
        % Fourier transform matrix
        F = exp( 2i*pi * u * (-fc:fc) );
        Gamma = compute_gamma(w,x,order);
        etaF = diag(conj(w)) * pinv(transpose(Gamma)) * [s; zeros(n*order,1)];
        eta = real( F * etaF );

        if 0
            %%% OLD %%%
            Gamma = PhiX(x,0);
            for k=1:order
                Gamma = [Gamma PhiX(x,k)];
            end
            p = pinv( Gamma )'*[s;zeros(n*order,1)];
            eta  = ConvolS(phi, p); % should be using bar phi
        end
        
        
    case 'fejer'
        
        % load phi operator for Candes            
        options.power = 4;
        phi = @(t,k)kernel('dirichlet', t, k, options);
        % evaluate the matrix Phi_x of size (P,N) at N=|x| sampling points
        Sampling = @(phi,x,k)phi( repmat(u, [1 n]) - repmat(x(:)', [P 1]), k );
        % convolves with reversed filter
        ConvolS = @(phi,x)real(ifft( conj(fft(phi(u,0))).*fft(x) ));
        % Phi_x matrices
        PhiX = @(x,k)Sampling(phi,x,k);
        % TODO: implement wih higher order derivatives ...
        X = repmat(x, [1 n]) - repmat(x', [n,1]);
        p = pinv([phi(X,0) phi(X,1); phi(X,1) phi(X,2)])*[s;zeros(n,1)];
        eta = [PhiX(x,0) PhiX(x,1)]*p;
        etaF = [];
        
    case 'positive'
        
        [X,U] = meshgrid(x, u);
        Z = sin(pi*(X-U)).^2;
        Z = prod( Z, 2);
        eta = 1 - 2*Z/max(Z);
        etaF = [];

    otherwise 
        error('Unknown certificate type');
end


function Gamma = compute_gamma(w,x,d)

% compute_gamma - compute Gamma_x matrix over Fourier domain
%
%   Gamma = compute_gamma(w,x,d);
%
%   Copyright (c) 2013 Gabriel Peyre

fc = (length(w)-1)/2;
x = x(:);
Ax = exp( 2i*pi * (-fc:fc)' * x' );

w1 = w(:);
Gamma = [];
for i=0:d
    Gamma = [Gamma, diag(w1) * Ax];
    % derivate the filter
    w1 = w1 .* 2i*pi .* (-fc:fc)';
end

    
