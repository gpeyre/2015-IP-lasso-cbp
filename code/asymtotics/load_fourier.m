function [Phi,Phi1,PhiD] = load_fourier(fc)

% load_fourier - load filtering operators over Fourier
%
%   [Phi,Phi1,PhiD] = load_fourier(fc);
%
% Phi is a partial Fourier matrix.
% Phi1 is the first derivative.
% PhiD is the kth derivative.
%
% Copyright (c) 2015 Gabriel Peyre

% operators
w = ones(2*fc+1,1);
Fourier = @(fc,x)exp(-2i*pi*(-fc:fc)'*x(:)');
Phi  = @(x)diag(w)*Fourier(fc,x);
% derivative
PhiD = @(x,k)diag(w)*diag( (-2i*pi*(-fc:fc)).^k )*Fourier(fc,x);
Phi1 = @(x)PhiD(x,1);