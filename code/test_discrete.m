%%
% Show the relationship between the solution on a discrete grid and on a
% continuous grid. 

save_output = 1;

%%
% Display options.

options.lw = 2;
options.msB = 30;
options.ms = 15;
options.ar = 1; % aspect ratio for plots

%%
% Parameters.

% number of measurements
fc = 6;
options.fc = fc;
% sampling point for display
P = 4096*4; 
options.P = P;
% Fourier coefficients of the filter
w = ones(2*fc+1,1);
KernelType = 'ideal';
% is the input measure on the grid ?
OnTheGrid = 1;
% dicretization step
Q = 128; 
% clear display
options.lgd = {};
Eta = {};

[Fourier,Phi,PhiS] = load_helpers();

%%
% Load the input measure.

% put it smaller to increase difficulty
if not(exist('delta'))
    delta = .8/fc;  % spacing between Diracs. 
end

DiracType = '3diracsc';
DiracType = '2diracsb';
DiracType = 'positive';
DiracType = 'evil';
DiracType = 'lots'; 
DiracType = '3diracsb';
DiracType = '3diracsa';
DiracType = '2diracsa';

rand('state', 123);
[x,s] = load_diracs(DiracType, delta, options);
a0 = s .* (1+rand(size(s)))/2; % value of the signal at the diracs

rep = ['results/discrete/' DiracType '/'];
if not(exist(rep))
    mkdir(rep);
end


%%
% Put the input measure on the grid or inbetween samples.

x = round(x*Q)/Q;
if not(OnTheGrid)
    % shift a bit to be off the grid
    x = x + (.2 + .6*rand(size(x)))/Q;
end
N = length(x);

%%
% Measurements.

y = Phi(w,x,a0);

%%
% Fuchs certificate on the Diracs' locations.

options.order = 0;
etaF0 = compute_certificate(w, x, s, options);
Eta{end+1} = etaF0;
options.lgd{end+1} = '\eta_F';

%%
% Check the location where Fuchs is not valid.

t = (find(abs(etaF0)>1)-1)/P;
[T,S] = meshgrid(t, x);
A = abs(T-S); A(A>1/2) = A(A>1/2) - 1;  % distance modulo 1
% if delta>1/Q, then the signal is identifiable in discrete, but not
% support-stable to noise
Dmin = max( min( abs(A) ) ); % max distance to nearest dirac
fprintf('If >1, not support stable: %3f\n', Q*Dmin);

%%
% Fuchs certificate on twin points

if OnTheGrid
    x1 = [x; x+1/Q];
else
    x1 = [floor(x*Q)/Q; ceil(x*Q)/Q];
end
s1 = [s;s];
% compute certificate
options.order = 0;
Eta{end+1} = compute_certificate(w, x1, s1, options);
options.lgd{end+1} = '\eta_0^{G}';

%%
% Vanishing certificate on the Diracs' locations.

options.order = 1;
Eta{end+1} = compute_certificate(w, x, s, options);
options.lgd{end+1} = '\eta_0';

%%
% Display.

options.ColorMode = 'default';
if save_output
    options.ColorMode = 'print';
    options.ColorMode = {'g-.' 'b--' 'r-'};
end

B = [0.4179 0.4492 0.9622 1.0372];  % box for zoom

clf;
% rectangle('Position', [B(1) B(3)  B(2)-B(1) B(4)-B(3)], 'FaceColor',[1 1 1]*.5);
plot_certificates(x, s, Eta, options);
drawnow;

if save_output
    str = [rep KernelType '-' DiracType '-' num2str(round(100*delta*fc)) '-ongrid' num2str(OnTheGrid)];
    saveas(gcf, str, 'epsc');
    fix_dottedline([str, '.eps']); % fix dashes in the .eps file
end

%%
% Zoom

axis(B);
% add the plot the discretized grid
plot( (0:Q-1)/Q, 1, 'm.', 'MarkerSize', 20 );
drawnow;
if save_output
    str = [rep KernelType '-' DiracType '-' num2str(round(100*delta*fc)) '-ongrid' num2str(OnTheGrid) '-zoom'];
    saveas(gcf, str, 'epsc');
    fix_dottedline([str, '.eps']); % fix dashes in the .eps file
end

return;

%%
% Solve the discrete Lasso problem.
% Vector should be identifiable.

% sampling matrix
Phi = Fourier( fc,(0:Q-1)/Q );
% primal solution
lambda = 0;
a = perform_lasso(y,Phi,lambda,1);
% primal solution
lambda = 1e-5;
a1 = perform_lasso(y,Phi,lambda,1);

%%
% Compute the regularization path.

Lambda = linspace(0,norm(Phi'*y, 'inf'),40);
A = [];
for i=1:length(Lambda)
    progressbar(i,length(Lambda));
    A(:,end+1) = perform_lasso(y,Phi,Lambda(i),1);
    clf;
    plot(log(abs(A(:,end))), '.-'); axis tight;
    drawnow;
end
I = find(sum(abs(A'))>1e-3);

clf;
plot(Lambda,A(I,:)', '-', 'LineWidth', 2);
axis tight; drawnow;
str = [rep KernelType '-' DiracType '-' num2str(round(100*delta*fc)) '-ongrid' num2str(OnTheGrid) '-path'];
saveas(gcf, str, 'epsc');


