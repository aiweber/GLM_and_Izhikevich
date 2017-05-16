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

cellType = 1; % type of Izhikevich neuron (numbered as in Izhikevich 2004)
%   Choose from:
%       1. tonic spiking
%       2. phasic spiking
%       3. tonic bursting
%       4. phasic bursting
%       5. mixed mode
%       6. spike frequency adaptation
%       7. Class 1
%       8. Class 2
%       9. spike latency
%       10. subthreshold oscillations -- not available
%       11. resonator
%       12. integrator
%       13. rebound spike
%       14. rebound burst
%       15. threshold variability
%       16. bistability
%       17. depolarizing after-potential -- not available
%       18. accomodation
%       19. inhibition-induced spiking
%       20. inhibition-induced bursting
%       21. bistability 2 (Not in original Izhikevich paper)

plotFlag = 1; % plot simulated response
saveFlag = 1; % save data to fid, in new folder
fid = pwd;    % root directory for project
T = 10000;    % max time (in ms)
jitter = 0;   % amount of jitter to add to spike times,
              %   uniformly distributed over [-jitter,jitter], measured in ms

[I, dt] = generate_izhikevich_stim(cellType,T);
[v, u, spikes, cid] = simulate_izhikevich(cellType,I,dt,jitter,plotFlag,saveFlag,fid);

%% STEP 2: fit GLM

% first choose parameters for basis vectors that characterize the
% stimulus and post-spike filters

%%% basis functions for stimulus filter
nkt = 100; % number of ms in stim filter
kbasprs.neye = 0; % number of "identity" basis vectors near time of spike;
kbasprs.ncos = 3; % number of raised-cosine vectors to use
kbasprs.kpeaks = [1 round(nkt/2)];  % position of first and last bump (relative to identity bumps)
kbasprs.b = 10; % how nonlinear to make spacings (larger -> more linear)
%%% basis functions for post-spike kernel
ihbasprs.ncols = 2;  % number of basis vectors for post-spike kernel
ihbasprs.hpeaks = [1 100];  % peak location for first and last vectors, in ms
ihbasprs.b = 100;  % how nonlinear to make spacings (larger -> more linear)
ihbasprs.absref = 0; % absolute refractory period, in ms

softRect = 0;    % use exponential nonlinearity
plotFlag = 1;    % plot fit
saveFlag = 1;    % save fit to fid, in new folder
maxIter = 1000;  % max number of iterations for fitting, also used for maximum number of function evaluations(MaxFunEvals)
tolFun = 1e-8;   % function tolerance for fitting
L2pen = 0;       % penalty on L2-norm of parameter coefficients

[k h dc prs kbasis hbasis] = fit_glm(I,spikes,dt,nkt,kbasprs,ihbasprs,[],softRect,plotFlag,maxIter,tolFun,L2pen);

% save
if saveFlag
    if ~isdir([fid '/glm_fits'])
        mkdir(fid, 'glm_fits')
    end
    tag = '';
    if softRect
        tag = '_sr';
    end
    save([fid '/glm_fits/' cid tag '_glmfit.mat'],'cellType','cid','dt','I','spikes','prs','kbasis','hbasis','softRect','maxIter','tolFun','L2pen','nkt','kbasprs','ihbasprs','k','h','dc')
    disp(['saved: ' fid '/glm_fits/' cid tag '_glmfit.mat'])
end


%% STEP 3: simulate responses of fit GLM

plotFlag = 1; % plot simulated data
saveFlag = 1; % save simulated data
runs = 10;    % number of trials to simulate

[y stimcurr hcurr r] = simulate_glm(I,dt,k,h,dc,runs,softRect,plotFlag);

if saveFlag
    if ~isdir([fid '/glm_sim_data'])
        mkdir(fid, 'glm_sim_data')
    end
    save([fid '/glm_sim_data/' cid tag '_glmdata.mat'],'y','stimcurr','hcurr','r')
    disp(['saved: ' fid '/glm_sim_data/' cid tag '_glmdata.mat'])
end


%% STEP 4: compare responses from simulated GLM and original Izhikevich data

compare_glm_to_iz(cellType,fid,softRect,jitter)

