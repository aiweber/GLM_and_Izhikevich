%% Example GLM fitting to Izhikevich neuron
%
% This code simulates data from an Izhikevich neuron, fits a GLM to it,
% simulates responses from that GLM, and then plots a comparison of the 
% simulated GLM responses to the original data.
%
% This code requires the following functions:
%   fminunc (in MATLAB's Optimization Toolbox)
%   generate_izhikevich_stim
%   simulate_izhikevich
%   fit_glm
%   simulate_glm
%   makeBasis_StimKernel
%   makeBasis_PostSpike
%   negloglike_glm_basis
%   negloglike_glm_basis_softrect
%   compare_glm_to_iz
%   normalizecols
%   sameconv
%   logexp1
% All except the first (fminunc) are provided in the package that contains
% this script.

%% STEP 1: simulate data from Izhikevich neuron

cellType = 3; % tonic bursting

plotFlag = 0; % plot simulated response
saveFlag = 1; % save data to fid, in new folder
fid = pwd;    % root directory for project
T = 3000;    % max time (in ms)
if cellType == 7 || cellType == 8
    T = 20000;  % these behaviors use multiple step heights, so generate more data
end
jitter = 0;   % amount of jitter to add to spike times,
              %   uniformly distributed over [-jitter,jitter], measured in ms

[I, dt] = generate_izhikevich_stim(cellType,T);



%%

Inew = gsmooth(randn(T/dt,1),1000)*200;
plot(Inew);

%%
[v1, u, spikes1, cid] = simulate_izhikevich(cellType,Inew,dt,jitter,plotFlag,saveFlag,fid);
plot(v1);
% [v1, u, spikes1, cid] = simulate_izhikevich(cellType,I*.5,dt,jitter,plotFlag,saveFlag,fid);
% [v2, u, spikes2, cid] = simulate_izhikevich(cellType,I*1,dt,jitter,plotFlag,saveFlag,fid);
% [v3, u, spikes3, cid] = simulate_izhikevich(cellType,I*2,dt,jitter,plotFlag,saveFlag,fid);
% tt = dt:dt:T;

%%
subplot(411);
plot(tt,I); set(gca, 'ylim', [0 12]);
subplot(412);
plot(tt, v1);
subplot(413);
plot(tt, v2);
subplot(414);
plot(tt, v3);

sps_train = [spikes1;spikes2;spikes3];
I_train = [I*0.5; I*1; I*2];


%% STEP 2: fit GLM


%%% basis functions for stimulus filter
nkt = 100; % number of ms in stim filter
kbasprs.neye = 0; % number of "identity" basis vectors near time of spike;
kbasprs.ncos = 7; % number of raised-cosine vectors to use
kbasprs.kpeaks = [.1 round(nkt/1.2)];  % position of first and last bump (relative to identity bumps)
kbasprs.b = 10; % how nonlinear to make spacings (larger -> more linear)
%%% basis functions for post-spike kernel
ihbasprs.ncols = 14;  % number of basis vectors for post-spike kernel
ihbasprs.hpeaks = [.1 100];  % peak location for first and last vectors, in ms
ihbasprs.b = 10;  % how nonlinear to make spacings (larger -> more linear)
ihbasprs.absref = 1; % absolute refractory period, in ms

softRect = 0;    % use exponential nonlinearity
plotFlag = 1;    % plot fit
saveFlag = 1;    % save fit to fid, in new folder
maxIter = 1000;  % max number of iterations for fitting, also used for maximum number of function evaluations(MaxFunEvals)
tolFun = 1e-12;  % function tolerance for fitting
L2pen = 0;       % penalty on L2-norm of parameter coefficients

[k, h, dc, prs, kbasis, hbasis] = fit_glm(I_train,sps_train,dt,nkt,kbasprs,ihbasprs,[],softRect,plotFlag,maxIter,tolFun,L2pen);

%% STEP 3: simulate responses of fit GLM

plotFlag = 0; % plot simulated data
saveFlag = 0; % save simulated data
runs = 10;    % number of trials to simulate

[y, stimcurr, hcurr, r] = simulate_glm(I_train,dt,k,h,dc,runs,softRect,plotFlag);


%%

clf;
tpl = 1:size(y,1);
plot(tpl, sps_train+runs+3);hold on;
for jj = 1:runs;
    plot(tpl,y(:,jj)+jj); 
end
hold off;


%% STEP 4: compare responses from simulated GLM and original Izhikevich data

compare_glm_to_iz(cellType,fid,softRect,jitter)

