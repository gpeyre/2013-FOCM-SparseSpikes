function [Fourier,Phi,PhiS] = load_helpers()

% load_helpers - load usefull functions
%   
% Fourier: transform on grid x for frequencies -fc:fc.
% Phi: compute y=Phi_{x}(a) where w are the weights. 
% PhiS: compute eta = Phi_u'(p), i.e. with result interpolated on grid u.
%
%   Copyright (c) 2013 Gabriel Peyre

Fourier = @(fc,x)exp(-2i*pi*(-fc:fc)'*x(:)');
Phi  = @(w,x,a)w .* (Fourier((length(w)-1)/2,x)*a);
PhiS = @(w,u,p)real( Fourier((length(w)-1)/2,u)'*(conj(w).*p) );