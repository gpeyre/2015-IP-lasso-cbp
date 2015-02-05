%%
% Test for the scaling of the constant involved in D-BP recovery.

addpath('toolbox_blasso/');
addpath('toolbox/');

rep = 'results/certificates/';
if not(exist(rep))
    mkdir(rep);
end

% cutoff frequency
fc = 10;

%% 
% Generate input measure.

% number of spikes
k = 2;
[x0,s0] = generate_spikes(k);
s0 = ones(k,1);

%%
% Display a few certificate for a given Delta.

delta_list = [1 .8 .6 .4 .2 .1];

for it=1:length(delta_list)
    % scale position
    delta = delta_list(it);
    x = .5 + delta*(x0-.5);    
    % compute certificates
    [etaV,etaW] = compute_certificate(x,s0,fc);
    % plot
    clf;
    plot_certificates(x, .8*s0, {etaV etaW});
    axis([0 1 -.4 1.1]);    
    saveas(gcf, [rep 'certificates-k' num2str(k) '-' num2str(it) '.eps'], 'epsc');
end
