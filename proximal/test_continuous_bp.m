%%
% Test for C-BP sparse spikes deconvolution, using a proximal FB algorithm.

addpath('toolbox/');

setfigname = @(name)set(gcf, 'Name', name, 'NumberTitle','off');

rep = 'results/';
if not(exist(rep))
    mkdir(rep);
end


%%
% Continuous kernel, defined on [0,1].

sigma = .01*3;
phi  = @(t)exp(-t.^2/(2*sigma^2));
phi1 = @(t)-t/(sigma^2).*exp(-t.^2/(2*sigma^2));

% Dirac's grid Grid size
N = 16*2;
Delta = 1/N;
% observation grid
rho = round(4096/N);
options.upsampling = rho;
P = N*rho;

%%
% Generate signals.

% Amplitude of displacement
if not(exist('dAmp'))
    dAmp = .95;
end
% Amplitude of elevations
aAmp = 1;

% number of spikes
k = 3;
% support
I = randperm(N); I = I(1:k);
% elevation
a0 = zeros(N,1);
a0(I) = aAmp * (1+rand(k,1))/2;
% displacement
d0 = zeros(N,1);
d0(I) = dAmp * 2*(rand(k,1)-.5);
%
b0 = d0.*a0;

% deterministic test
I = round( N/(2*k):N/k:N );
a0 = zeros(N,1); a0(I) = [.6 1 .8];
d0 = zeros(N,1); d0(I) = [-.2 1 -.7] * dAmp;
b0 = d0.*a0;

%%
% Generate observations.

% callbacks
[Gamma,GammaS,PhiExact,Phi,PhiS,Psi,PsiS] = load_filters(phi,phi1,N,options);
% exact observation
x0 = (0:N-1)'/N + d0*Delta/2;
y = PhiExact(x0,a0);
% basis pursuit approximation (quantification)
y0 = Phi(a0);
% C-BP approximation (linear approximation)
y1 = Gamma([a0 b0]);

figure(1); clf;
subplot(3,1,1);
plot(y); axis tight; title('Observations');
subplot(3,1,2);
plot([y-y0 y-y1]); axis tight; title('0th order');
legend('0th order', '1st order');
subplot(3,1,3);
plot(y-y1, 'g'); axis tight; title('1st order');


%%
% Solve C-BP.

% regularization parameter
lambda = 1e-5;
lambda = 60;

% Observations
Y = y;
if 0
    % this should return exact values
    Y = y1;
    lambda = 1e-8;
end

options.niter = 4000;
options.fbdamping = 1.8;
options.method = 'fista';
options.method = 'fb';
[a,b,delta,x, R] = perform_cbp(phi,phi1,Y,lambda,N,options);


% aspect ratio for plots
fs = 25; % font size
ar = 9/16; % aspect ratio
lw = 2; % line width for plots
msB = 30;
setDisp = @()set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', fs);
mystem = @(x,y, col)stem(x, y, [col '.--'], 'MarkerSize', msB, 'LineWidth', lw);
% set(gca, 'XTick', [], 'YTick', [-1 1]);

%%
% Display.

t = (0:N-1)'/N;
s = (0:P-1)'/P;
J = find(a>1e-3);

figure(1); setfigname('Recovery');
clf; hold on;
plot(s, y, 'LineWidth', 2);
mystem(x0(I), a0(I), 'k'); % initial spikes
mystem(t(J) + delta(J), a(J), 'r');  % recovered spikes
axis([0 1 0 1]);
setDisp(); box on;

figure(2); setfigname('Recovery (a,b)'); clf;
subplot(2,1,1);
hold on;
stem(a0, 'k.', 'MarkerSize', 20);
plot(a, 'r.'); axis tight;
subplot(2,1,2);
hold on;
stem(b0, 'k.', 'MarkerSize', 20);
plot(b, 'r.'); axis tight;

if 0
    figure(2); setfigname('Energy decay'); clf;
    sel = 1:options.niter/4;
    plot(sel, log(R(sel)/min(R)-1), '.-'); axis tight;
end

%%
% Compute a full homotopy.

% compute the homotopy path
Lmax = max(PhiS(y));
Llist = linspace(Lmax,1e-6,30);

% warmup stage
options.initialization = zeros(N,2);
options.niter = 5000;
[a,b,delta,x, R] = perform_cbp(phi,phi1,Y,Llist(1),N,options);
%
options.niter = 500;
options.initialization = [a,b];
A = []; B = [];
for i=1:length(Llist)
    [a,b,delta,x, R] = perform_cbp(phi,phi1,Y,Llist(i),N,options);
    A(:,end+1) = a; B(:,end+1) = b;
    options.initialization = [a,b];
end


Ic= setdiff(1:N,I);
figure(3); setfigname('Homotopy'); clf;
subplot(2,1,1); hold on;
plot(Llist, A(I,:)', 'LineWidth', lw);
plot(Llist, A(Ic,:)', 'k:', 'LineWidth', lw); axis tight; box on;
subplot(2,1,2); hold on;
plot(Llist, B(I,:)', 'LineWidth', lw);
plot(Llist, B(Ic,:)', 'k:', 'LineWidth', lw); axis tight; box on;
axis([0 max(Llist) -1 1]);
%
v = round(100*max(abs(d0)));
str = [rep 'homotopy-' num2str(v)];
saveas(gcf, str, 'epsc');
fix_dottedline([str, '.eps']); % fix dashes in the .eps file
