function [v, u, spikes, cid] = simulate_izhikevich(cellType,I,dt,jitter,plotFlag,saveFlag,fid)
% [v u spikes cid] = simulate_izhikevich(cellType,I,dt,jitter,plotFlag,saveFlag,fid)
%
%  This code generates data from an Izhikevich neuron of a type specified
%  by the user and saves it to a .mat file in the specified directory.
%
%  The inputs are:
%   cellType: type of Izhikevich neuron.  Choose from:
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
%   I: stimulus (current input)
%   dt: time step, in ms  
%   jitter: add jitter to spike times, uniformly distributed over [-jitter,jitter], 
%           measured in ms 
%           (Will not actually change output of Izhikevich model, but will jitter 
%           the timing of spikes in output vector 'spikes'.)
%   plotFlag: 0 or 1, plot simulated data
%   saveFlag: 0 or 1, save simulated data to file
%   fid: root directory for project, defaults to current directory
%   
%  The outputs are:
%     v: voltage response of the neuron
%     u: membrane recovery variable
%     spikes: vector of 0s and 1s indicating spikes
%     cid: string that identifies cell type

%% set defaults

if ~exist('jitter','var') || isempty(jitter)
    jitter = 0;
end

if ~exist('plotFlag','var') || isempty(plotFlag)
    plotFlag = 1;
end

if ~exist('saveFlag','var') || isempty(saveFlag)
    saveFlag = 1;
end

if ~exist('fid','var') || isempty(fid)
    fid = pwd;
end


%% parameters
% numbered as in Izhikevich 2004
% a-d values taken from code published by Eugene Izhikevich, found at:
% http://www.izhikevich.org/publications/izhikevich.m
% only 1-9, 11-16, 18-21 are suitable for GLM


%        a         b       c      d        
pars = [0.02      0.2     -65     6;      % 1. tonic spiking
        0.02      0.25    -65     6;      % 2. phasic spiking
        0.02      0.2     -50     2;      % 3. tonic bursting
        0.02      0.25    -55     0.05;   % 4. phasic bursting
        0.02      0.2     -55     4;      % 5. mixed mode
        0.01      0.2     -65     5;      % 6. spike frequency adaptation
        0.02      -0.1    -55     6;      % 7. Class 1
        0.2       0.26    -65     0;      % 8. Class 2
        0.02      0.2     -65     6;      % 9. spike latency
        0.05      0.26    -60     0;      % 10. subthreshold oscillations
        0.1       0.26    -60     -1;     % 11. resonator
        0.02      -0.1    -55     6;      % 12. integrator
        0.03      0.25    -60     4;      % 13. rebound spike
        0.03      0.25    -52     0;      % 14. rebound burst
        0.03      0.25    -60     4;      % 15. threshold variability
        1         1.5     -60     0;      % 16. bistability
        1         0.2     -60     -21;    % 17. depolarizing after-potential
        0.02      1       -55     4;      % 18. accomodation
        -0.02     -1      -60     8;      % 19. inhibition-induced spiking
        -0.026    -1      -45     0;      % 20. inhibition-induced bursting
        1         1.5     -60     0 ];    % 21. bistability 2 (Not in original Izhikevich paper)

cids = {'RS' 'PS' 'TB' 'PB' 'MM' 'FA' 'E1' 'E2' 'SL' 'SO' 'R' 'I' 'ES' 'EB' 'TV' 'B' 'DA' 'A' 'IS' 'IB' 'B2'};

%%%
a = pars(cellType,1);
b = pars(cellType,2);
c = pars(cellType,3);
d = pars(cellType,4);
cid = cids{cellType};

T = length(I)*dt;
t = dt:dt:T;


%% initialize variables

threshold = 30;

v = zeros(length(t),1);
u = zeros(length(t),1);
spikes = zeros(length(t),1);

% different initial v and u values to start different neuron types near
% stable fixed point (prevent spiking in absence of inputs near t=0)
if sum([16 21] == cellType) % if bistable
    v(1) = -54;
    u(1) = -77; 
elseif cellType == 12 % integrator
    v(1) = -90;
    u(1) = 0;
elseif sum([19 20]==cellType) % inhibition-induced spiking/bursting
    v(1) = -100;
    u(1) = 80;
else
    v(1) = -70; 
    u(1) = -14; 
end

% Izhikevich model doesn't show this kind of bistability, so simulate
% responses using first form of bistability
Iplot = I;
if cellType == 21
    I = abs(I+65)-65;
end


%% run model

for tt = 1:length(I)-1
    dvdt = 0.04*v(tt)^2 + 5*v(tt) +140 - u(tt) + I(tt);
    v(tt+1) = v(tt) + dvdt*dt;
    
    dudt = a*(b*v(tt+1)-u(tt));
    u(tt+1) = u(tt) + dudt*dt;
    
    if v(tt+1)>threshold
        v(tt) = threshold;  % makes spikes of uniform height
        v(tt+1) = c;
        u(tt+1) = u(tt+1) + d;
        spikes(tt+1) = 1;
    end
end



%% if jitter ~= 0, add noise to spike times
if jitter  % add noise to spike times
    spikeIdx = find(spikes);
    jitters = round((rand(size(spikeIdx))-.5)*2 * jitter/dt);
    spikeIdx = spikeIdx + jitters;
    spikes = zeros(size(spikes));
    spikes(spikeIdx) = 1;
end


%% plot results

if plotFlag
    I = Iplot;
    sTimes = t(logical(spikes));
    
    warning('off','MATLAB:hg:ColorSpec_None')
    figure;
    
    subplot(2,1,1)
    plot(t,I)
    ylabel('current')
    ylim([min(I)-.05*abs(min(I)) max(I)+.05*abs(max(I))])
    xlim([dt T])
    box off
    title('stimulus')
    
    subplot(2,1,2); hold on;
    plot(t,v)
    for s = 1:length(sTimes)
        plot([sTimes(s) sTimes(s)],[threshold*1.05 threshold*1.05+.2*threshold],'k')
    end
    xlim([dt T])
    xlabel('time (ms)'); ylabel('voltage')
    title('response')
    
end

%% save stimulus/response
if saveFlag
    if ~isdir([fid '/izhikevich_data'])
        mkdir(fid, 'izhikevich_data')
    end
    save([fid '/izhikevich_data/' cid '_iz.mat'],'cellType','cid','dt','T','a','b','c','d','threshold','I','u','v','spikes')
    disp(['saved: ' fid '/izhikevich_data/' cid '_iz.mat'])
end