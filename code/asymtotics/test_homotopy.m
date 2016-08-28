%%
% Test for C-BP resolution using CVX.
% Change cvx_path variable in install_cvx if needed.

addpath('toolbox/');
addpath('toolbox_blasso/');
install_cvx();

rep = 'results/homotopy/';
if not(exist(rep))
    mkdir(rep);
end

%%
% Helpers.

fs = 20; % font size
set_tick = @()set(gca, 'XTick', 0:10:100, 'YTick', -2:.5:2, 'FontSize', fs);
set_ar = @()set(gca, 'PlotBoxAspectRatio', [1 1/3 1]);
setfigname = @(name)set(gcf, 'Name', name, 'NumberTitle','off');

% select the problem to be tested.
method = 'bp';
method = 'cbp';

% sampling for display purpose.
P = 2048;
t = (0:P-1)'/P;
% operators
fc = 11;
Fourier = @(fc,x)exp(-2i*pi*(-fc:fc)'*x(:)');
Phi  = @(x)Fourier(fc,x);
% derivative
PhiD = @(x,k)diag( (-2i*pi*(-fc:fc)).^k )*Fourier(fc,x);
Phi1 = @(x)PhiD(x,1);

%%
% Generate input measure.

% grid for solving
N = 256*2;
N = 256;
Delta = 1/N;
z = (0:N-1)'/N;

% number of spikes
k = 2;
[x0,s0,q] = generate_spikes(k);

switch k
    case 1
        delta = 1; 
    case 2
        delta = .22; % validate C-BP source condition
        delta = .22;
        delta = .2; % invalidate C-BP source condition
    case 3
        delta = .65;
    otherwise
        error('You need to setup k');
end
x0 = (x0-.5)*delta + .5; % do the scaling
s0 = [1 .6 .75]; s0 = s0(1:k);

% data exactly on the grid
x0 = round(x0*N)/N;
I = round(x0*N); % position on grid
a0 = zeros(N,1); a0(I) = s0;

% observations (here noiseless)
y = Phi(z) * a0;

str0 = ['homotopy-k' num2str(k) '-d' num2str(round(delta*100))];
str = [str0 '-n' num2str(N) '-' method];

% vizualizes the filtered signal
clf; hold on;
u = real(Phi(t)'*y); u = u/max(u);
plot(t, u);
stem(x0, s0, 'k');
axis tight;

%%
% Plot certificate

[etaV,etaW] = compute_certificate(x0,s0,fc);
clf;
plot_certificates(x0, .8*sign(s0), {etaV etaW});
axis([0 1 -.4 1.1]);
saveas(gcf, [rep str0 '-certificates.eps'], 'epsc');

%%
% Solve C-BP or BP. The recovered spikes have locations
%       z + Delta/2*b/a

% derivative matrix
switch method 
    case 'bp'
        G = zeros(size(Phi(z)));
    case 'cbp'
        G = Delta/2*Phi1(z);
end

lambda = 0;
options.verb = 'quiet';
options.precision = 'best';
[a,b] = solve_cbp(Phi(z), G,y,lambda, options);

q = 600;
q = 100*2;
zoom = 1;
switch N
    case 128
        lambda_max = 40;
    case 256
        lambda_max = 10; % no zoom
        if zoom==1
            lambda_max = .004; % super zoom for C-BP
        end
    case 512
        lambda_max = 20; % full path
        lambda_max = 10; % zoom on start    
        if zoom==1
            lambda_max = .01/2; % super zoom for C-BP
        end
    case 1024
        lambda_max = 2; % zoom on start
    otherwise
        error('Need to setup lambda_max');
end
lambda_list = linspace(0,lambda_max,q);
A = []; B = [];
for i=1:q
    progressbar(i,q);
    lambda = lambda_list(i);
    [A(:,i),B(:,i)] = solve_cbp(Phi(z), G,y,lambda, options);
end


J = find( sum(A')>1e-3 ); % display only active components
lw = 2;

figure(1); setfigname('A vectors');
clf;
plot_homotopy(lambda_list, A(J,:), J, I, z(J));
set_ar(); set_tick(); axis tight; box on;
saveas(gcf, [rep str '-A-z' num2str(zoom) '.eps'], 'epsc');

if strcmp(method, 'cbp')
    
figure(2); setfigname('B vectors');
clf;
plot_homotopy(lambda_list, B(J,:), J, I, z(J));
set_ar(); set_tick(); axis tight; box on;
saveas(gcf, [rep str '-B-z' num2str(zoom) '.eps'], 'epsc');

figure(3); setfigname('B/A vectors');
clf;
C = B(J,:)./A(J,:);
i = find(A(J,:)<1e-3);
C(i) = NaN;
% disambiguates
t = linspace(.8, 1, size(C,1))';
C = C .* repmat(t,[1 q]);
plot_homotopy(lambda_list, C, J, I, z(J));
set_ar(); set_tick(); axis tight; box on;
saveas(gcf, [rep str '-BA-z' num2str(zoom) '.eps'], 'epsc');

end
