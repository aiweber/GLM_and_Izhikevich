function [y stimcurr hcurr r] = simulate_glm(x,dt,k,h,dc,runs,softRect,plotFlag)
% [y stimcurr hcurr r] = simulate_glm(x,dt,k,h,dc,runs,softRect,plotFlag)
%
%  This code fits a Poisson GLM to given data, using basis vectors to
%  characterize the stimulus and post-spike filters.
%
%  The inputs are:
%   x: stimulus
%   dt: time step of x and y in ms
%   k: stimulus filter
%   h: post-spike filter
%   dc: dc offset
%   runs: number of trials to simulate
%   softRect: 0 uses exponential nonlinearity; 1 uses soft-rectifying nonlinearity
%   plotFlag: 0 or 1, plot simulated data
%
%  The outputs are:
%   y: spike train (0s and 1s)
%   stimcurr: output of stimulus filter (without DC current added)
%   hcurr: output of post-spike filter
%   r: firing rate (stimcurr + hcurr + dc passed through nonlinearity)
%
%  This code requires the function sameconv, as well as logexp1 if using soft-rectified nonlinearity.

%% set defaults

if ~exist('runs','var') || isempty(runs)
    runs = 5;
end

if ~exist('softRect','var') || isempty(softRect)
    softRect = 0;
end

if ~exist('plotFlag','var') || isempty(plotFlag)
    plotFlag = 0;
end

%% generate data with fitted GLM

nTimePts = length(x);
refreshRate = 1000/dt; % stimulus in ms, sampled at dt

if softRect
    NL = @logexp1;
else
    NL = @exp;
end

g = zeros(nTimePts+length(h),runs); % filtered stimulus + dc
y = zeros(nTimePts,runs); % initialize response vector (pad with zeros in order to convolve with post-spike filter)
r = zeros(nTimePts+length(h)-1,runs); % firing rate (output of nonlinearity)
hcurr = zeros(size(g));

stimcurr = sameconv(x,k);
Iinj = stimcurr + dc;

for runNum = 1:runs
    
    g(:,runNum) = [Iinj; zeros(length(h),1)]; % injected current includes DC drive
    
    %%% loop to get responses, incorporate post-spike filter
    for t = 1:nTimePts
        r(t,runNum) = feval(NL,g(t,runNum));  % firing rate
        if rand<(1-exp(-r(t,runNum)/refreshRate)); % 1-P(0 spikes)
            y(t,runNum) = 1;
            g(t:t+length(h)-1,runNum) = g(t:t+length(h)-1,runNum) + h;  % add post-spike filter
            hcurr(t:t+length(h)-1,runNum) = hcurr(t:t+length(h)-1,runNum) + h;
        end
    end
end

hcurr = hcurr(1:nTimePts,:);  % trim zero padding
r = r(1:nTimePts,:);  % trim zero padding

%
if plotFlag
    
    minT = 1/dt;
    maxT = length(x);
    
    tIdx = minT:maxT;
    t = (tIdx-minT)*dt;
    
    %%% stimulus
    subplot(4,1,1); hold on;
    plot(t,x(tIdx),'linewidth',2)
    xlim([min(t) max(t)])
    ylim([min(x(tIdx))-.05*abs(min(x(tIdx))) max(x(tIdx))+.05*abs(max(x(tIdx)))])
    box off
    title('stimulus')
    
    %%% filter outputs
    subplot(4,1,2); hold on;
    plot(t,Iinj(tIdx),'linewidth',1.5);
    plot(t,hcurr(tIdx,1),'r','linewidth',1.5)
    xlim([min(t) max(t)])
    ylim([min([hcurr(tIdx,1); Iinj(tIdx)])*1.1 max([hcurr(tIdx,1); Iinj(tIdx)])*1.1])
    box off
    title('filter outputs')
    
    %%% prob(firing)
    subplot(4,1,3); hold on;
    semilogy(t,feval(NL,hcurr(tIdx,1)+Iinj(tIdx)),'color',[.5 .5 .5],'linewidth',1.5)
    xlim([min(t) max(t)])
    box off
    title('prob(firing)')
    
    %%% GLM spikes
    subplot(4,1,4); hold on;
    spikeHeight = .7;
    
    for i = 1:size(y,2) % for each run of glm simulation
        spt = find(y(tIdx,i));
        for spikeNum = 1:length(spt)
            plot([spt(spikeNum)*dt spt(spikeNum)*dt],[i-.5 i-.5+spikeHeight],'color',[.5 .5 .5],'linewidth',1.25)
        end
    end
    
    xlim([0 max(t)-min(t)])
    ylim([0 runs+spikeHeight])
    xlabel('time (ms)')
    title('spikes')
    
end
