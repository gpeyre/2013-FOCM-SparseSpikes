function [x,s] = load_diracs(DiracType, delta, options)

% load_diracs - load sum of diracs from a name
%
%   [x,s] = load_diracs(DiracType, delta, options);
%
%   delta is the spacing between Diracs
%   x are the positions of the Diracs
%   s are the corresponding signs
%
%   Copyright (c) 2013 Gabriel Peyre

options.null = 0;
fc = getoptions(options, 'fc', -1);

switch DiracType
    case '2diracsa'
        x = [.5-delta/2 .5+delta/2]';
        s = [1 -1]';
    case '2diracsb'
        x = [.5-delta/2 .5+delta/2]';
        s = [1 1]';
    case '2diracsMod'
        % 2 close Diracs and 2 far away
        x = [.15 .5-delta/2 .5+delta/2 .85]';
        s = [1 1 -1 1]';
    case '3diracsa'
        x = [.5-delta .5 .5+delta]';
        s = [1 1 -1]';
    case '3diracsb'
        x = [.5-delta .5 .5+delta]';
        s = [1 -1 1]';
    case '3diracsc'
        x = [.5-delta .5 .5+delta]';
        s = [1 1 1]';
    case 'evil'
        delta = 1.01/fc;
        s = [+1,+1,+1,-1,+1,+1,+1]';
        n = length(s);
        x = delta*(0:n-1)'; x = x-mean(x)+.5;
    case 'lots'
        % regular sampling
        if fc<0
            error('You need to provide options.fc');
        end            
        delta = 1.1/fc;
        % rng(11); % for fc=26 delta=1.01  -> 5
        x = (0:delta:1-1/fc)'; 
        n = length(x);
        s = sign(randn(n, 1));
    case 'lotsb'
        % iregular sampling
        if fc<0
            error('You need to provide options.fc');
        end
        delta = 1.1/fc;
        % rng(11); % for fc=26 delta=1.01  -> 5
        a = rescale(rand(1000,1).^2,1,5)/fc;
        x = [0; cumsum(a)];
        x = x(x<1-1/fc); 
        n = length(x);
        s = sign(randn(n, 1));
    case 'tvtest'
        % An example for TV
        s = .13; % No staircasing for fc=6
        s = .3;
        x = [1/6 1/2 1/2+s]';
        s = [1 1 -1]';
    case 'positive'
        k = getoptions(options, 'nb_diracs', 3); % number of diracs
        x0 = 1/2; % centeral points
        a = linspace(-1/2,1/2,k)'; % uniform
        % delta = randn(k,1);
        x = 1/2 + delta*a*k;
        s = ones(k,1);
        %% The following ones are for fc=26
    case 'pathological1'
        x = [0.0034,  0.0525,  0.1011,  0.1501,  0.1992,  0.2514,  0.3033,  0.3555,  0.406 ,  0.4541,  0.5051,  0.5533,  0.6049,  0.6562,  0.7052,  0.754 ,  0.8032,  0.8525,  0.9036,  0.9519];
        s = [ 1., -1.,  1., -1., -1.,  1., -1.,  1.,  1., -1., -1., -1., -1.,  1.,  1., -1., -1.,  1., -1., -1.];
    case 'pathological2'
        x = [ 0.0082,  0.0566,  0.1097,  0.1651,  0.2208,  0.2717,  0.3206,  0.3701,  0.4266,  0.4803,  0.5358,  0.592 ,  0.6408,  0.6895,  0.7431,  0.7975,  0.8495,  0.9022,  0.9519];
        s = [ 1.,  1., -1.,  1.,  1.,  1., -1., -1.,  1., -1.,  1., -1.,  1., -1.,  1., -1., -1.,  1., -1.];
    case 'pathological3'
        x = [ 0.0019,  0.0631,  0.1126,  0.1714,  0.2207,  0.2711,  0.3218,  0.3779,  0.4401,  0.4908,  0.5401,  0.6006,  0.658 ,  0.7207,  0.7802,  0.8304,  0.8878,  0.9519];
        s = [-1.,  1., -1.,  1., -1., -1.,  1., -1.,  1., -1.,  1., -1.,  1., -1.,  1.,  1., -1.,  1.];
        x = 1-x;
    otherwise
        error('Unknown');
end
s = sign(s(:)); x = x(:);