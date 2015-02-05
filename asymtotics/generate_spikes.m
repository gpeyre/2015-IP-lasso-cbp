function [x0,s0,q] = generate_spikes(k)

% generate_spikes - generate a bunch of k spikes
%
%   [x0,s0,q] = generate_spikes(k); 
%
% Copyright (c) 2014 Gabriel Peyre

switch k
    case 1;
        x0 = .5; 
        s0 = 1; 
        q = 1; 
    case 2
        x0 = [.3 .7]';
        s0 = [1 1]';
        % extended support direction
        q = [-1 +1]'; % direction
    case 3
        x0 = [.25 .4 .65]';
        s0 = [1 1 1]';
        % extended support direction
        q = [-1 +1 -1]'; % direction
    case 4
        x0 = [.23 .43 .55 .77]';
        x0 = (1/(2*k):1/k:1)';
        x0 = x0 + [+1 -1 +1 -2]' * 1/(4*k);
        s0 = [1 1 1 1]';
        % extended support direction
        q = [+1 -1 +1 -1]'; % direction
    otherwise
        x0 = (1/(2*k):1/k:1)' + 1/(3*k) * randn(k,1);
        s0 = ones(k,1);
        q = ones(k,1);
end