%%
% Display regularization paths (positions/values as a function of lambda).

addpath('toolbox/');

rep = 'results/paths/';
if not(exist(rep))
    mkdir(rep);
end

%%
% Helpers 

[Fourier,Phi,PhiS] = load_helpers();

%%
% Display options.


% display options
options.lw = 2;
options.msB = 30;
options.ms = 15;
options.ar = 1; % aspect ratio for plots
options.fs = 25; % font size

% cuttoff frequency
if not(exist('fc'))
    fc = 26;
    fc = 10;
end
options.fc = fc;


w = ones(2*fc+1,1);
options.w = w;


%%
% Load input measure.

if not(exist('DiracType'))
    DiracType = 'lotsb';
    DiracType = 'pathological2';
    DiracType = '2diracsa';
    DiracType = '3diracsa';
end

% spacing between Diracs. 
switch DiracType
    case '2diracsa'
        delta = .6/fc;  
    case '3diracsa'
        delta = .7/fc; 
    otherwise
        delta = .9/fc; 
end

[x,s] = load_diracs(DiracType, delta, options);

% value of the signal at the diracs
a0 = s;
if length(x)==3
    a0 = [.6 1 -1]';
elseif length(x)>3
    a0 = s .* (1+rand(size(s)))/2;
end

%%
% Display dual certificates. 

options.order = 1;
options.type = 'leastsquare';
etaV = compute_certificate(w, x, s, options); 
Eta = etaV;
options.lgd = {'\eta_V'};

options.ColorMode = 'default';
clf;
plot_certificates(x, a0, Eta, options);
drawnow;
saveas(gcf, [rep DiracType '-fc' num2str(fc) '-certif.eps'], 'epsc');

%%
% Input noise.

if not(exist('sigma'))
    sigma = .5; % for 2 diracs
end
randn('state', 123);
Noise = (randn(2*fc+1,1)+1i*randn(2*fc+1,1))*sigma;

%%
% Observations.

y = Fourier(fc,x)*a0 + Noise;
% observation over spacial domain
u = Fourier(fc,linspace(0,1,2048))'*y;

%%
% Compute solutions for varying lambda.

K = 1500; % sampling density
lambda_list = linspace(0,max(abs(u))/2, K);
x1 = {}; a1 = {};
for i=1:length(lambda_list)
    progressbar(i,length(lambda_list));
    lambda = lambda_list(i);
    [x1{i},a1{i}] = perform_sdp_superresolution(y,lambda);
end

%%
% Display the solution path. 

smax = 30; % size max of dots
clf; hold on;
for i=1:round(length(x1))
    for j=1:length(x1{i})
        v = a1{i}(j);
        % width
        m = abs(v)^.5 * smax;
        % color
        c = [max(v,0), 0, max(-v,0)];
        c = min(5*c,1);
        % plot
        plot( lambda_list(i), x1{i}(j), 'k.', 'MarkerSize', m, 'MarkerFaceColor', c, 'MarkerEdgeColor', c );
    end
end
box on; axis tight;
xlabel('\lambda');
axis([0 max(lambda_list) 0 1]);
set(gca, 'YTick', [0 .5 1], 'FontSize', options.fs);
saveas(gcf, [rep DiracType '-fc' num2str(fc) '-noise' num2str(round(10*sigma)) '-path.eps'], 'epsc');
