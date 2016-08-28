%%
% Test for resolution of the dual problem using Douglas-Rachford splitting.

addpath('./toolbox/');

% number of sampling points
P = 2048*2;
options.P = P;
u = (0:P-1)'/P;
% display options
options.lw = 2;
options.msB = 30;
options.ms = 15;
options.ar = 1; % aspect ratio for plots
save_output = 1; % save as .eps files the output (no legend)

% cuttoff frequency
if not(exist('fc'))
    fc = 6; 
end
options.fc = fc;

%% 
% Load the kernel

KernelType = 'ideal'; % ideal low pass one
w = ones(2*fc+1,1);
options.w = w;


%% 
% Signal parameters.

% put it smaller to increase difficulty
if not(exist('delta'))
    delta = .7/fc;  % spacing between Diracs. 
end

DiracType = '3diracsb';
DiracType = '3diracsc';
DiracType = '2diracsb';
DiracType = 'evil';
DiracType = 'lotsb';
DiracType = 'lots';
DiracType = '3diracsa';

rep = ['results/certificates/' KernelType '/' DiracType '/'];
if not(exist(rep))
    mkdir(rep);
end

[x,s] = load_diracs(DiracType, delta, options);


Fourier = @(fc,x)exp(-2i*pi*(-fc:fc)'*x(:)');
y = w .* (Fourier(fc,x)*s);
lambda =  .05; % 1e-3;


%%
% Minimal norm certificate, computed using CVX.

options.solver = 'cvx';
options.h = w;
p0 = perform_sdp_dual(y,lambda, options);
eta0 = real( Fourier(fc,u)'*(conj(w).*p0) ); % High-res Fourier matrix

%%
% Solve using DR.

options.solver = 'dr';
p1 = perform_sdp_dual(y,lambda, options);
eta1 = real( Fourier(fc,u)'*(conj(w).*p1) ); % High-res Fourier matrix

clf;
plot([eta0, eta1]);
axis tight;
legend('IP', 'DR');