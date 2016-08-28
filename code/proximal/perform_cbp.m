function [a,b,delta,x, R] = perform_cbp(phi,phi1,y,lambda,N,options)

% perform_cbp - solve the continuous basis-pursuit deconvolution method
%
%   Solve the problem
%       min_{a,b} 1/2*|y - Phi(a) - Delta/2 * Phi'(b)|^2 + lambda*|a|_1
%   s.t.   a>=0  and  |b|<=a.
%
%   y is the observed signal of dimension P=rho*N. 
%   Delta=1/N is the grid spacing (uniform) on the Dirac's grid.
%   a is the amplitude of the estimated spikes.
%   delta=Delta/2 * b/a is the displacement
%   x = (i/N + delta_i)_i are the recovered spikes location
%
%   rho=options.upsampling is an integer telling how much the observed
%   signal is upsampled (rho=1 implies the same grid is used for Dirac's
%   location and observations). 
%
%   phi and phi1 are callback returning the value of phi and phi' at any point in [0,1] 
%   (should be able to take a vector of t as input).
%
%   The resolution is performed using the forward-backward algorithm.
%
% This method is detailed in:
%   Recovery of Sparse Translation-Invariant Signals With Continuous Basis Pursuit
%   Ekanadham, C. ; Tranchina, D. ; Simoncelli, E.P.
%   IEEE Transactions on Signal Processing, 59(10), 2011
%
%   You can set in options the parameter of perform_fb
%
%   Copyright (c) 2014 Gabriel Peyre

options.null = 0;
rho = getoptions(options, 'upsampling', 1);
P = N*rho;
Delta = 1/N;

% initialization 
u0 = getoptions(options, 'initialization', zeros(N,2));

[Gamma,GammaS,PhiExact,Phi,PhiS,Psi,PsiS] = load_filters(phi,phi1,N,options);

%%
% Projection on C. 

C  = @(w)w(:,2)+1i*w(:,1);
Ci = @(v)[imag(v), real(v)];
ProjOct = @(v)max(real(v),0) + 1i*max(imag(v),0);
ProjC = @(w)Ci(exp(1i*pi/4)*ProjOct(exp(-1i*pi/4)*C(w)));

%% 
% Prox of J(w)=lambda*|x|_1+i_C(x,s)

ProxJ = @(w,lambda)ProjC( w-[lambda*ones(size(w,1),1) zeros(size(w,1),1)] );

%% 
% Solve for
%    min_w 1/2*|y-A(w)|^2 + lambda*J(w)

% estimate Lipschitz constant of nablaF
[L,e] = compute_operator_norm(@(w)GammaS(Gamma(w)),randn([N 2]));

F = @(w)1/2*norm(y-Gamma(w), 'fro')^2;
gradF = @(w)GammaS(Gamma(w)-y);
E = @(w)F(w) + lambda*norm(w(:,1),1);
options.report = @(w)E(w); 
[u,R] = perform_fb(u0, @(w,tau)ProxJ(w,tau*lambda), gradF, L, options);
%
a = u(:,1); b = u(:,2);
delta = Delta/2 * b./a; delta(a<1e-9) = 0;
x = (0:N-1)'/N + delta;

end
