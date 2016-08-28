%%
% Evaluation of the different pre-certificates.

addpath('./toolbox/');

% save or not in .eps file the output
save_output = 0; 

%%
% Helpers 

[Fourier,Phi,PhiS] = load_helpers();

%% 
% Parameters.

% number of sampling points for the display of the pre-certificates.
P = 2048*8;
options.P = P;
u = (0:P-1)'/P;
% display options
options.lw = 2;
options.msB = 30;
options.ms = 15;
options.ar = 1; % aspect ratio for plots

% cuttoff frequency
if not(exist('fc'))
    fc = 6; 
end
options.fc = fc;

%% 
% Load the kernel

KernelType = 'ideal';
w = ones(2*fc+1,1);
options.w = w;

%% 
% Signal parameters.

% put it smaller to increase difficulty
if not(exist('delta'))
    delta = .7/fc;  % spacing between Diracs. 
end

% select here the configuration you want to bench.
DiracType = '3diracsb';
DiracType = '3diracsc';
DiracType = '2diracsb';
DiracType = 'evil';
DiracType = 'lotsb';
DiracType = 'lots';
DiracType = 'pathological1'; % eta_0 is saturating outside support, eta_V'' has wrong sign
DiracType = 'pathological2'; % eta_V'' has wrong sign
DiracType = 'pathological3';
DiracType = '3diracsa';

rep = ['results/certificates/' DiracType '/'];
if not(exist(rep))
    mkdir(rep);
end

[x,s] = load_diracs(DiracType, delta, options);
Eta = {}; options.lgd = {};

%%
% Measurements (noiseless).

y = Phi(w,x,s);

%%
% Fuchs pre-certif.

if not(save_output)
    options.order = 0;
    options.type = 'leastsquare';
    Eta{end+1} = compute_certificate(w, x, s, options);
    options.lgd{end+1} = '\eta_F';
end


%%
% Minimal norm certificate.

lambda =  0.005;
options.solver = 'dr'; % use proximal method
options.niter = 500; % only useful for DR
options.gamma = 1/100; % only useful for DR
options.solver = 'cvx'; % use interor point method
p = perform_sdp_dual(y,lambda, options);
eta0 = PhiS(w,u,p);
Eta{end+1} = eta0;
options.lgd{end+1} = '\eta_0';

%%
% Candes/Fernandez-Granda pre-certificate.

options.order = 1;
options.type = 'fejer';
Eta{end+1} = compute_certificate(w, x, s, options);
options.lgd{end+1} = '\eta_{CF}';

%%
% Vanishing derivative pre-certif.

options.order = 1;
options.type = 'leastsquare';
etaV = compute_certificate(w, x, s, options); 
Eta{end+1} = etaV;
options.lgd{end+1} = '\eta_V';

%%
% Plot the different certificates.

options.ColorMode = 'default';
if save_output
    options.ColorMode = 'print';
    options.ColorMode = {'g-.' 'b--' 'r-'};
end
clf;
plot_certificates(x, s, Eta, options);
drawnow;

%%
% Save results.

if save_output
    str = [rep KernelType '-' DiracType '-' num2str(round(100*delta*fc))];
    saveas(gcf, str, 'epsc');
    fix_dottedline([str, '.eps']); % fix dashes in the .eps file
end



w = 0.2771-0.1299;
h = 2.2;
clf;
rectangle('Position', [0.122 -h/2 w h], 'FaceColor',[1 1 1]*.5);
plot_certificates(x, s, Eta, options);
saveas(gcf, str, 'epsc');
fix_dottedline([str, '.eps']); % fix dashes in the .eps file

clf;
options.lgd = '';
plot_certificates(x, s, Eta, options);
axis([0.1299    0.2771   -1.1429    1.1221]);
saveas(gcf, [str '-zoom'], 'epsc');
fix_dottedline([str '-zoom.eps']); % fix dashes in the .eps file

