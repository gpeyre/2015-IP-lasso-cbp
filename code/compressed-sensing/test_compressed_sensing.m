%%
% Test for the size of the extended support in a compressed sensing scenario.
% We compute the minimal norm certificate using DR (Douglas-Rachford).

addpath('toolbox');

rep = 'results/';
[~,~,~]=mkdir(rep);

%%
% Useful helpers

mynorm = @(x)norm(x(:));
% set aspect ratio
set_ar = @()set(gca, 'PlotBoxAspectRatio', [1 1/3 1]);
% set figure title
setfigname = @(name)set(gcf, 'Name', name, 'NumberTitle','off');

% dimension
n = 400;
% Number of measurements.
p = round(n/4);
% Measurement matrix
A = randn(p,n)/sqrt(p);

% proximal map of L1
proxL1 = @(x,gamma)max(0,1-gamma./max(1e-15,abs(x))).*x;
% projector on y=Phi*x
pA = A'*(A*A')^(-1);
projPrimal = @(x,y)x + pA*(y-A*x);

%%
% Primal problem:   min F(x)+G(x) = |x|_1 + iota_{C}(x)   C={x ; A*x=A*x0}

% number of trials
q = 4000;
% List of benched sparsity.
slist = 14:2:42; % Ok for perfect recovery
slist = 1:2:42;
[x0,Slist] = gen_sparse_signals(slist,n,q);

% Measurements.
y = A*x0;

% Perform DR on the set of signals |x0|. Note that the proximal mappings
% operate in parallel on all the signals in |x0|.
% Each |i|, count the average number |proba(i)|
% of recovered vectors of sparsity |slist(i)| (up to a given, small, precision).
lun = [];
options.niter = 2000;
options.mu = 1;
options.gamma = 1;
options.verb = 1;
options.report = @(x)mynorm(A*x-y); % should be 0
% monitor the L1 norm
options.report = @(x)norm(x(:), 1);
% run DR algorithm
[x,R] = perform_dr(x0,proxL1,@(x,tau)projPrimal(x,y),options);

% Compute IC criteria
C = mat2cell(x0,n,ones(q,1));
IC = cellfun(@(x)compute_ic_l1(x,A),C);

% compute proba of exact recovery
probaID = [];
probaIC = [];
EID = mean(abs(x-x0))<.05;
EIC = IC<1;
for j=1:length(slist)
    s = slist(j);
    probaID(j) = mean(EID(Slist==s));
    probaIC(j) = mean(EIC(Slist==s));
end

% display proba of recovery
lw = 2;
clf; hold on;
plot(slist, probaIC, 'k--', 'LineWidth', lw);
plot(slist, probaID, 'k', 'LineWidth', lw);
axis([min(slist) max(slist) -.03 1.03]);
box on;
saveas(gcf, [rep 'ident-vs-ic.eps'], 'epsc');


%%
% Minimal solution to dual problem:   min F(p,eta)+G(p,eta) = 1/2*|p|^2 + iota_{C}(eta) + iota_D(p,eta)
%   C={eta ; eta_I=sign(x_0I), |eta_{I^c}|_inf<=1}
%   D={p,eta ; eta=A^* p}
%     proj_D(p,eta) = ( eta0, A^* eta0 ), eta0=U*(A^* eta + p), U=(A^*A + Id)^{-1}

% Targeted sparsity
slist = [10 12 14 16 18];

for s=slist

    q = 4000;
    [x0,Slist] = gen_sparse_signals(s,n,q);

    % extractor
    P = @(v)v(1:p,:);
    Eta = @(v)v(p+1:n+p,:);

    % impose sign and L^inf
    ProjLinf = @(eta)clamp(eta,-1,1);
    ProxS = @(eta) eta .* (abs(x0)<1e-10) +  sign(x0);
    ProjLinfS = @(eta)ProjLinf(ProxS(eta));

    % 1/2*|p|^2 + iota_{C}(eta)
    ProxF = @(v,gamma)[ P(v)/(1+gamma); ProjLinfS(Eta(v)) ];
    % iota_D(p,eta), proj_D(p,eta) = ( p0, A^* p0 ), p0=U*(A^* eta + p)
    U = (A*A' + eye(p))^(-1);
    ProxG = @(v,gamma)[ U*(A*Eta(v) + P(v)) ; A'*U*(A*Eta(v) + P(v)) ];

    % run DR algorithm
    options.report = @(x)0;
    [v,R] = perform_dr(zeros(n+p,q),ProxF,ProxG,options);
    eta = Eta(v);

    % extract saturating points
    tol = .01;
    Sat = abs(abs(eta)-1)<tol; % saturation points
    % number of saturations
    Ext = sum(Sat);
    % plot histogram
    elist = s-1:max(Ext);
    elist = 12:51;
    h = hist(Ext,elist);
    % plot
    clf;
    hold on;
    bar(elist(1:end-1), h(1:end-1), 'k');
    plot([s s], [0 max(h(1:end-1))], 'r:', 'LineWidth', 2);
    axis tight;
    set(gca, 'YTick', []);
    set_ar();
    box on;
    saveas(gcf, [rep 'distrib-extsup-s' num2str(s) '.eps'], 'epsc');

end
