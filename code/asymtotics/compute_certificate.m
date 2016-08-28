function [etaV,etaW] = compute_certificate(x,s0,fc)

% compute_certificate - compute BLASSO and C-BP certificates
%
% [etaV,etaW] = compute_certificate(x,s0,fc); 
%
% etaV is for BLASSO
% etaW is for C-BP
%
% Copyright (c) 2015 Gabriel Peyre

% sampling for display purpose.
P = 2048*4;
t = (0:P-1)'/P;

[Phi,Phi1,PhiD] = load_fourier(fc);
s0 = sign(s0(:));
k = length(s0);

% Vanishing derivative
pV = pinv( [Phi(x) PhiD(x,1)])' * [s0; zeros(k,1)];
etaV = real( Phi(t)' * pV );
% Vanishing 3rd derivative
pV3 = pinv( [Phi(x) PhiD(x,1) PhiD(x,3)])' * [s0; zeros(2*k,1)];
etaW = real( Phi(t)' * pV3 );