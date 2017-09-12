%%
% This code recreates Figure 7 of Weber & Pillow 2017 ("Capturing the
% dynamical repertoire of single neurons with generalized linear models")
%
% It sweeps over a range of two different parameter values of the GLM.  One
% controls the amplitude of a single basis vector of the stimulus filter;
% the other controls the amplitude of a single basis vector of the
% post-spike filter.
%
% After simulating responses from each GLM, the code automatically
% classifies each response as quiescent, phasic spiking, phasic bursting,
% tonic spiking, or tonic bursting. The classification results, along with
% example filters, are then plotted.
%
% This code requires the following functions:
%   makeBasis_StimKernel
%   makeBasis_PostSpike
%   sameconv
%   fitdist
%   chi2gof
%   gmdistribution
% The last 3 are part of MATLAB's Statistics Toolbox.  All others are
% provided in the package that contains this script.

%%
clearvars
saveDir = pwd; % directory to save files

%% make step stimulus

x = zeros(140000,1);
x(2000:120000) = .4;
dt = .1;

%% create basis functions and initialize parameters

nkt = 100;  % 100, number of ms in stim filter

%%% arbitrary number of basis functions for stimulus filter
kbasprs.neye = 0; % Number of "identity" basis vectors near time of spike;
kbasprs.ncos = 2; % Number of raised-cosine vectors to use
kbasprs.kpeaks = [0 round(nkt/4)];  % Position of first and last bump (relative to identity bumps)
kbasprs.b = 100;                    % Offset for nonlinear scaling (larger -> more linear)
kbasisTemp = makeBasis_StimKernel(kbasprs,nkt);
nkb = size(kbasisTemp,2);
lenkb = size(kbasisTemp,1);
kbasis = zeros(lenkb/dt,nkb);
for bNum = 1:nkb
    kbasis(:,bNum) = interp1((1:lenkb)',kbasisTemp(:,bNum),linspace(1,lenkb,lenkb/dt)');
end

%%% basis functions for post-spike kernel
ihbasprs.ncols = 2;         % Number of basis vectors for post-spike kernel
ihbasprs.hpeaks = [10 40];  % Peak location for first and last vectors
ihbasprs.b = 50;            % How nonlinear to make spacings
ihbasprs.absref = 0;        % Absolute refractory period
[ht,hbas,hbasis] = makeBasis_PostSpike(ihbasprs,dt);
hbasis = [zeros(1,ihbasprs.ncols); hbasis]; % enforce causality: post-spike filter only affects future time points

nkbasis = size(kbasis,2); % number of basis functions for k
nhbasis = size(hbasis,2); % number of basis functions for h

prs = zeros(nkbasis+nhbasis+1,1); % initialize parameters

xconvki = zeros(size(x,1),nkbasis);
yconvhi = zeros(size(x,1),nhbasis);

for knum = 1:nkbasis
    xconvki(:,knum) = sameconv(x,kbasis(:,knum));
end


%% simulate glm
prsStims = -1.75:.25:.25;
prsHists = -1:.25:1;
numReps = 25;
spAll = cell(numReps,length(prsStims),length(prsHists));
ks = cell(length(prsStims),1);
hs = cell(length(prsHists),1);
for repNum = 1:numReps
    
    stimK = 0;
    for prsStim = prsStims(1:length(prsStims))
        stimK = stimK + 1;
        disp(prsStim)
        
        histK = 0;
        for prsHist = prsHists(1:length(prsHists))
            histK = histK + 1;
            
            %%%% parameters
            prs = [1 prsStim prsHist -3 -1]';
            
            NL = @exp;
            
            nTimePts = length(x);
            nkbasis = size(kbasis,2);
            refreshRate = 1000/dt; % stimulus in ms, sampled at dt from izhikevich model
            
            k = kbasis*prs(1:nkbasis); % k basis functions weighted by given parameters
            h = [-1000*ones(round(5/dt),1); hbasis*prs(nkbasis+1:end-1)]; h(1) = 0; % h basis functions weighted by given parameters, with refractory period
            dc = prs(end); % dc current (accounts for mean spike rate)
            
            if repNum == 1 && stimK == 1
                hs{histK} = h;
            end
            if repNum == 1 && histK == 1
                ks{stimK} = k;
            end
                
            g = zeros(nTimePts+length(h),1); % filtered stimulus + dc
            y = zeros(nTimePts,1); % initialize response vector (pad with zeros in order to convolve with post-spike filter)
            r = zeros(nTimePts+length(h)-1,1); % firing rate (output of nonlinearity)
            hcurr = zeros(size(g));
            
            xconvk = sameconv(x,k);
            Iinj = xconvk + dc;
            
            g = [Iinj; zeros(length(h),1)]; % injected current includes DC drive
            
            %%% loop to get responses, incorporate post-spike filter
            for t = 1:nTimePts
                r(t) = feval(NL,g(t));  % firing rate (generator signal through nonlinearity)
                if rand<(1-exp(-r(t)/refreshRate))
                    y(t) = 1;
                    g(t:t+length(h)-1) = g(t:t+length(h)-1) + h;  % add post-spike filter
                    hcurr(t:t+length(h)-1) = hcurr(t:t+length(h)-1) + h;
                end
            end
            sp = find(y);
            spAll{repNum,stimK,histK} = sp;
                        
        end
    end
    
end
disp('DONE WITH GLM SIMULATION')

save([saveDir '/glm_param_sweep_example'],'spAll','nkt','kbasprs','ihbasprs','kbasis','hbasis','dc','ks','hs','x','dt','saveDir')

%% classify behavior

minIdx = find(x,1,'first');
maxIdx = find(x,1,'last');
classesAll = zeros(length(prsStims),length(prsHists),numReps);

for repNum = 1:numReps
    
    stimK = 0;
    for prsStim = prsStims(1:length(prsStims))
        stimK = stimK + 1;
        
        histK = 0;
        for prsHist = prsHists(1:length(prsHists))
            histK = histK + 1;
            
            sp = spAll{repNum,stimK,histK};
            y = zeros(size(x));
            y(sp) = 1;
            
            ISIs = diff(sp(sp>=minIdx & sp<=maxIdx));
            uISIs = unique(ISIs);
            
            if length(uISIs)>30  % need to have a large number of unique ISIs
                pd = fitdist(ISIs,'exponential');
                chi2 = chi2gof(ISIs,'CDF',pd,'alpha',.01);
            else 
                chi2 = 1; 
            end
            
            if sum(y(2000:4000))~=0 && sum(y(20000:120000))<5 % spikes at short latency, very few at longer latency
                if sum(y(2000:10000))< 3
                    classi = 2;  % phasic spiking
                else
                    classi = 5; % phasic bursting
                end
                
            elseif sum(y(2000:4000))==0 && sum(y(20000:120000))<5 % few spikes
                classi = 1;  % quiescent
                
            elseif length(unique(ISIs)) == 1 % perfectly regular spiking
                classi = 3;  % regular spiking
                
            else  % compare AIC for uni- and bi-modality
                
                f1 = gmdistribution.fit(ISIs,1,'regularize',1e-6);
                f2 = gmdistribution.fit(ISIs,2,'regularize',1e-6);
                
                if f1.AIC*.9 < f2.AIC % unimodal
                    classi = 3;  % regular spiking
                else
                    classi = 4;  % bursting
                end
            end
            
            classesAll(stimK,histK,repNum) = classi;
            
        end
    end
end

disp('DONE WITH CLASSIFICATION')
save([saveDir '/glm_param_sweep_example'],'classesAll','spAll','nkt','kbasprs','ihbasprs','kbasis','hbasis','dc','ks','hs','x','dt','saveDir')

%% plot results

figure;
subplot(5,5,[2:5 7:10 12:15 17:20]); hold on;
r = [120 65 64 72 85 99 127 181 217 230 230 217]'/255;
g = [28 59 101 139 161 173 185 189 173 142 100 33]'/255;
b = [129 147 177 194 177 153 185 76 60 52 44 32]'/255;
baseMap = [r g b];
sDiffs = unique(diff(prsStims));
hDiffs = unique(diff(prsHists));
for i = 2:length(prsStims)
    for j = 1:length(prsHists)
        prsStim = prsStims(i);
        prsHist = prsHists(j);
        
        rectangle('position',[prsHist-.5*hDiffs prsStim-.5*sDiffs hDiffs sDiffs],'facecolor',baseMap(mode(classesAll(i,j,:)),:),'edgecolor',baseMap(mode(classesAll(i,j,:)),:));
    end
end
xlim([min(prsHists)-hDiffs*.51 max(prsHists)+hDiffs*.51]);
ylim([min(prsStims)+sDiffs*.51 max(prsStims)+sDiffs*.51]);
set(gca,'xtick',prsHists(1:2:end),'ytick',prsStims(2:2:end))
xlabel({'post-spike filter'; 'component 1 amplitude'})
ylabel({'stimulus filter'; 'component 2 amplitude'})

subplot(5,5,21); hold on;
labels = {'quiescent','phasic spiking','tonic spiking','tonic bursting','phasic bursting'};
for i = 1:5
    text(0,i/3,labels{i},'color',baseMap(i,:))
end
ylim([.2 2])
axis off

subplot(5,5,1); hold on;
plot(ks{end},'linewidth',2,'color','b')
plot([0 length(ks{end})],[0 0],'--','color',[1 1 1]*.6);
ylim([-.25 .35])

subplot(5,5,16); hold on;
plot(ks{1},'linewidth',2,'color','b')
plot([0 length(ks{1})],[0 0],'--','color',[1 1 1]*.6);
ylim([-.25 .35])

subplot(5,5,22); hold on;
plot(hs{1},'linewidth',2,'color','r')
plot([0 length(hs{1})],[0 0],'--','color',[1 1 1]*.6);
ylim([-5 0.1])
xlim([0 length(hs{1})])

subplot(5,5,25); hold on;
plot(hs{end},'linewidth',2,'color','r')
plot([0 length(hs{end})],[0 0],'--','color',[1 1 1]*.6);
ylim([-5 0.1])
xlim([0 length(hs{end})])
