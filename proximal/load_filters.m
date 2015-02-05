function [Gamma,GammaS,PhiExact,Phi,PhiS,Psi,PsiS] = load_filters(phi,phi1,N,options)

% load_filters - load various operators
%
%   [Gamma,GammaS,PhiExact,Phi,PhiS,Psi,PsiS] = load_filters(phi,phi1,N,options);
%
%   phi and phi1 are callbacks which implement phi(t) and phi'(t) for
%   arbitrary t in [0,1]. 
%
%   Copyright (c) 2014 Gabriel Peyre

rho = getoptions(options, 'upsampling', 1);
P = N*rho;
Delta = 1/N;

% convolution and adjoint of convolution
convol  = @(x,h)real(ifft(fft(x).*fft(h)));
convolS = @(x,h)real(ifft(fft(x).*conj(fft(h))));

% filters on the upsampled grid.
t = [0:P/2, -P/2+1:-1]' / P;
phi_d = phi(t);
phi1_d = phi1(t);
% operators
Phi  = @(x)convol(upsample(x,rho),phi_d);
PhiS = @(x)downsample(convolS(x,phi_d),rho);
Psi  = @(s)convol(upsample(s,rho),phi1_d)*Delta/2;
PsiS = @(s)downsample(convolS(s,phi1_d),rho)*Delta/2;
Gamma  = @(w)Phi(w(:,1)) - Psi(w(:,2));
GammaS = @(y)[PhiS(y), -PsiS(y)];
certify_adjoint(Gamma,GammaS,[N,2],1e-9);

PhiExact = @(x0,alpha0)SumSpikes(phi,x0,alpha0,t);

end

%%
function y = SumSpikes(phi,x0,alpha0,t)

% true signal
y = zeros(length(t),1);
for i=1:length(x0)
    T = t-x0(i); T = mod(T,1); T(T>.5) = T(T>.5)-1;
    y = y + alpha0(i) * phi( T );
end

end

%%
function y = upsample(x,rho)

y = zeros(length(x)*rho,1);
y(1:rho:end) = x;


end

function x = downsample(y,rho)

x = y(1:rho:end);

end