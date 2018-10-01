function [k h dc prs kbasis hbasis] = fit_glm(x,y,dt,nkt,kbasprs,ihbasprs,prs,softRect,plotFlag,maxIter,tolFun,L2pen)
% [k h dc prs kbasis hbasis] = fit_glm(x,y,dt,nkt,kbasprs,ihbasprs,prs,softRect,plotFlag,maxIter,tolFun,L2pen)
%
%  This code fits a Poisson GLM to given data, using basis vectors to
%  characterize the stimulus and post-spike filters.
%
%  The inputs are:
%   x: stimulus
%   y: spiking data, vector of 0s and 1s
%   dt: time step of x and y in ms
%   nkt: number of ms in stimulus filter
%   kbasprs: structure containing parameters of stimulus filter basis vectors
%       kbasprs.neye: number of "identity" basis vectors near time of spike
%       kbasprs.ncos: number of raised-cosine vectors to use
%       kbasprs.kpeaks: position of first and last bump (relative to identity bumps)
%       kbasprs.b: how nonlinear to make spacings (larger -> more linear)
%   ihbasprs: structure containing parameters of post-spike filter basis vectors
%       ihbasprs.ncols: number of basis vectors for post-spike kernel
%       ihbasprs.hpeaks: peak location for first and last vectors
%       ihbasprs.b: how nonlinear to make spacings (larger -> more linear)
%       ihbasprs.absref: absolute refractory period, in ms
%   prs: vector to initialize fit parameters
%   softRect: 0 uses exponential nonlinearity; 1 uses soft-rectifying nonlinearity
%   plotFlag: 0 or 1, plot simulated data
%   maxIter: maximum number of iterations for fitting
%   tolFun: function tolerance for fitting
%   L2pen: size of L2 penalty on coefficients in prs (defaults to 0)
%
%  The outputs are:
%   k: stimulus filter
%   h: post-spike filter
%   dc: DC offset
%   prs: full set of coefficients for basis vectors, [k_coeffs h_coeffs dc]
%   kbasis: basis vectors for stimulus filter
%   hbasis: basis vectors for post-spike filters
%
%  This code requires the function fminunc in the MATLAB Optimization
%  Toolbox, as well as the following functions:
%   makeBasis_StimKernel
%   makeBasis_PostSpike
%   normalizecols
%   sameconv
%   negloglike_glm_basis (or negloglike_glm_basis_softrect, if using soft-rectified nonlinearity)
%   logexp1 (if using soft-rectified nonlinearity)

%% set defaults

if ~exist('nkt','var') || isempty(nkt)
    nkt = 100;
end

if ~exist('kbasprs','var') || isempty(kbasprs)
    %%% basis functions for stimulus filter
    kbasprs.neye = 0; % number of "identity" basis vectors near time of spike;
    kbasprs.ncos = 3; % number of raised-cosine vectors to use
    kbasprs.kpeaks = [1 round(nkt/2)];  % position of first and last bump (relative to identity bumps)
    kbasprs.b = 10; % how nonlinear to make spacings (larger -> more linear)
end

if ~exist('ihbasprs','var') || isempty(ihbasprs)
    %%% basis functions for post-spike kernel
    ihbasprs.ncols = 2;  % number of basis vectors for post-spike kernel
    ihbasprs.hpeaks = [1 100];  % peak location for first and last vectors,in ms
    ihbasprs.b = 10;  % how nonlinear to make spacings (larger -> more linear)
    ihbasprs.absref = 0; % absolute refractory period, in ms
end

if ~exist('softRect','var') || isempty(softRect)
    softRect = 0;
end

if ~exist('plotFlag','var') || isempty(plotFlag)
    plotFlag = 0;
end

if ~exist('maxIter','var') || isempty(maxIter)
    maxIter = 100;
end

if ~exist('tolFun','var') || isempty(tolFun)
    tolFun = 1e-8;
end

if ~exist('L2pen','var') || isempty(L2pen)
    L2pen = 0;  % penalty on L2 norm
end

refreshRate = 1000/dt; % stimulus in ms, sampled at dt


%% create basis functions and initialize parameters

kbasisTemp = makeBasis_StimKernel(kbasprs,nkt);
nkb = size(kbasisTemp,2);
lenkb = size(kbasisTemp,1);
kbasis = zeros(lenkb/dt,nkb);
for bNum = 1:nkb
    kbasis(:,bNum) = interp1([1:lenkb]',kbasisTemp(:,bNum),linspace(1,lenkb,lenkb/dt)');
end

[ht,hbas,hbasis] = makeBasis_PostSpike(ihbasprs,dt);
hbasis = [zeros(1,ihbasprs.ncols); hbasis]; % enforce causality: post-spike filter only affects future time points

nkbasis = size(kbasis,2); % number of basis functions for k
nhbasis = size(hbasis,2); % number of basis functions for h

if ~exist('prs','var') || isempty(prs)
    prs = zeros(nkbasis+nhbasis+1,1); % initialize parameters
end

%%

xconvki = zeros(size(y,1),nkbasis);
yconvhi = zeros(size(y,1),nhbasis);

for knum = 1:nkbasis
    xconvki(:,knum) = sameconv(x,kbasis(:,knum));
end

for hnum = 1:nhbasis
    yconvhi(:,hnum) = sameconv(y,flipud(hbasis(:,hnum)));
end

%% minimization

warning('off','optim:fminunc:SwitchingMethod')
opts = optimoptions('fminunc','algorithm','trust-region','gradobj','on','hessian','on','display','iter','maxiter',maxIter,'maxfunevals',maxIter,'tolfun',tolFun,'tolX',tolFun);
if softRect
    NL = @logexp1;
    fneglogli = @(prs) negloglike_glm_basis_softRect(prs,NL,xconvki,yconvhi,y,1,refreshRate);
else
    NL = @exp;
    fneglogli = @(prs) negloglike_glm_basis(prs,NL,xconvki,yconvhi,y,1,refreshRate,L2pen);
end

%% optimization
prs = fminunc(fneglogli,prs,opts);

%% calculate filters from basis fcns/weights
k = kbasis*prs(1:nkbasis); % k basis functions weighted by given parameters
h = hbasis*prs(nkbasis+1:end-1); % k basis functions weighted by given parameters
dc = prs(end); % dc current (accounts for mean spike rate)

%% plot results
if plotFlag
    figure;
    subplot(2,2,1); hold on;
    for i = 1:size(kbasis,2)
        plot(kbasis(:,i))
    end
    xlim([1 length(k)])
    set(gca,'xtick',0:round(length(k)/5):length(k),'xticklabel',round(fliplr(-dt*(0:round(length(k)/5):length(k)))))
    
    subplot(2,2,3)
    plot(k)
    xlim([1 length(k)])
    set(gca,'xtick',0:round(length(k)/5):length(k),'xticklabel',round(fliplr(-dt*(0:round(length(k)/5):length(k)))))
    xlabel('time (ms)')
    title('stimulus filter')
    
    subplot(2,2,2); hold on;
    for i = 1:size(hbasis,2)
        plot(hbasis(:,i))
    end
    xlim([1 length(h)])
    set(gca,'xtick',0:round(length(h)/5):length(h),'xticklabel',round(dt*(0:round(length(h)/5):length(h))))
    subplot(2,2,4)
    plot(h)
    xlim([1 length(h)])
    set(gca,'xtick',0:round(length(h)/5):length(h),'xticklabel',round(dt*(0:round(length(h)/5):length(h))))
    xlabel('time (ms)')
    title('post-spike filter')
end

