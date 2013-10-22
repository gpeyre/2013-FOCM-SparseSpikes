function plot_certificates(x, s, Eta, options)

% plot_certificates - display of Diracs and certificates
%
%   plot_certificates(x, s, Eta, options);
%
%   x is position of Diracs
%   s is signs of Diracs
%   Eta{i} is a certificate
%
%   options.lgd{i} is legend for certificate #i
%
%   Copyright (c) 2013 Gabriel Peyre

if not(iscell(Eta))
    Eta = {Eta};
end


options.null = 0;
lgd = getoptions(options, 'lgd', {});
msB = getoptions(options, 'msB', 30);
ms = getoptions(options, 'ms', 15);
lw = getoptions(options, 'lw', 2);
ar = getoptions(options, 'ar', 1/3);
fs = getoptions(options, 'fs', 25);
ar = [1 ar 1];

%%
% Set color for display of certificates.

ColorMode = getoptions(options, 'ColorMode', 'default');
q = length(Eta);
if iscell(ColorMode)
    Col = ColorMode;
    ColorMode = '';
end
Mkr = {};
switch ColorMode
    case 'default'
        Col = {'r' 'b' 'g' 'm' 'c'};
        if q>length(Col)
            error('Set options.ColorMode to jet mode.');
        end
    case 'print'
        Col = { 'k--' 'b' 'r' 'm' 'c--'};
        Mkr = {'none', '.', 'none', '.'};
        ms = 15;
        if q>length(Col)
            error('Set options.ColorMode to jet mode.');
        end
    case 'jet'
        C = jet(q);
        for i=1:q
            Col{i} = C(i,:);
        end
    case ''
        % do nothing
    otherwise 
        error('Unknown color mode');
end
        
P = length(Eta{1});
u = (0:P-1)'/P;

hold on;
% plot support point
stem(x, s, 'k.--', 'MarkerSize', msB, 'LineWidth', 2);
plot([0 1],  [1 1], 'k--', 'LineWidth', 2); 
plot([0 1], -[1 1], 'k--', 'LineWidth', 2);
axis([0 1 -1.2 1.2]);
set(gca, 'XTick', [], 'YTick', [-1 1], 'PlotBoxAspectRatio', ar, 'FontSize', fs);
box on;
% plot the pre-certificate
h = [];
for i=1:length(Eta)
    h(end+1) = plot(u, Eta{i}, Col{i}, 'LineWidth', lw);
end
% add markers
if not(isempty(Mkr))
    Q = 60; % number of markers
    t = round( linspace(1,length(u), Q) );
    for i=1:length(Eta)
        plot(u(t), Eta{i}(t), Col{i}, 'Marker', Mkr{i}, 'MarkerSize', ms, 'LineStyle', 'none');
    end
end

if not(isempty(lgd))
    k = min([length(h) length(lgd)]);
    legend(h(1:k), {lgd{1:k}});
end