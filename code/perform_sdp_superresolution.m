function [x1,a1] = perform_sdp_superresolution(y, lambda)

% perform_sdp_superresolution - perform super-resolution over Radon measures
%
%   [x1,a1] = perform_sdp_superresolution(y, lambda, w);
%
%   The operator is
%       Phi(m) = m * phi   where   hat phi = 1_{[-fc,fc]}
%   The reconstruction for lambda=0 is
%       min_{Phi(m)=y} |m|_TV    where   y=Phi(m_{a,x})
%       where  m_{a,x} = sum_k a_k delta_{x_k}
%   The reconstruction for lambda>0 is
%       min_m  1/2*|Phi(m)-y|^2 + lambda*|m|_TV    where   y=Phi(m_{a,x})
%
%   This code is based on the code provided in
%       Towards a Mathematical Theory of Super-Resolution
%       Emmanuel J. Candes and Carlos Fernandez-Granda
%
%   Copyright (c) 2013 Gabriel Peyre

n = length(y);
fc = (n-1)/2;

if nargin<2
    lambda = 0;
end

if 0
x = x(:)';
a = a(:);
% data 
F = exp(-2i*pi*k'*x); % Fourier matrix
y = F*a + w; 
end


%% Solve SDP
u = perform_sdp_dual(y,lambda);

%% Roots of dual polynomial
aux_u =- conv(u,flipud(conj(u)));
aux_u(2*fc+1)=1+aux_u(2*fc+1);
roots_pol = roots(flipud(aux_u));

if 0
    clf;
    plot(real(roots_pol),imag(roots_pol),'*');
    hold on;
    plot(cos(2*pi*x), sin(2*pi*x),'ro');
    plot( exp(1i*linspace(0,2*pi,200)), '--' );
    hold off;
    legend('Roots','Support of x'); axis equal; axis([-1 1 -1 1]*1.3);
end

% Isolate roots on the unit circle
tol = 1e-2;
roots_detected = roots_pol(abs(1-abs(roots_pol)) < tol);

% BUG
% [auxsort,ind_sort]=sort(real(roots_detected));
% BUG
[auxsort,ind_sort]=sort(angle(roots_detected));
roots_detected = roots_detected(ind_sort);
roots_detected = roots_detected(1:2:end);

% Roots are double so take 1 and out of 2 and compute argument
x1 = angle(roots_detected)/2/pi;
% Argument is between -1/2 and 1/2 so convert angles
x1 = sort(mod(x1,1));
% Accuracy 
% delta = abs(sort(x1') - sort(x));

k = -fc:fc;
F_est = exp(-2i*pi*k'*x1'); % estimated Fourier matrix
s0 = sign(real(F_est'*u));
a1 = real(F_est\y - lambda*pinv(F_est'*F_est)*s0 ); % estimated amplitudes
