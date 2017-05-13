function compare_glm_to_iz(cellType,fid,softRect,jitter)
% compare_glm_to_iz(cellType,fid,softRect,jitter)
%
%  This function loads simulated responses of an Izhikevich neuron and a
%  GLM fit to that neuron and plots them for comparison.
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
%   fid: root directory for project, defaults to pwd
%   softRect: 0 uses exponential nonlinearity; 1 uses soft-rectifying nonlinearity
%   jitter: amount of jitter added to spike times, measured in ms 
%
%  This code requires the function logexp1 if using soft-rectified nonlinearity.

%% set defaults

if ~exist('softRect','var') || isempty(softRect)
    softRect = 0;
end

if ~exist('jitter','var') || isempty(jitter)
    jitter = 0;
end


%% plot GLM simulation results against Izhikevich data

Ts = [400 1100;
    450 1250;
    900 1300;
    950 1250;
    950 1250;
    350 1250;
    9250 15250;
    8250 14250;
    750 1050;
    200 1100;
    1900 4000;
    5100 9100;
    600 900;
    600 900;
    1400 2200;
    125 225;
    200 1100;
    300 2200;
    200 1100;
    200 1100;
    525 725;
    100 1000;
    100 10000];


%% load simulated data
cids = {'RS' 'PS' 'TB' 'PB' 'MM' 'FA' 'E1' 'E2' 'SL' 'SO' 'R' 'I' 'ES' 'EB' 'TV' 'B' 'DA' 'A' 'IS' 'IB' 'B2' 'RS2' 'SP'};
cid = cids{cellType};
if jitter
    cid = [cid '_jitter'];
end
tag = '';
if softRect
    tag = [tag '_sr'];
    NL = @logexp1;
else
    NL = @exp;
end

load([fid '/glm_fits/' cid tag '_glmfit.mat'])
load([fid '/glm_sim_data/' cid tag '_glmdata.mat'])
load([fid '/izhikevich_data/' cid '_iz.mat'])

%%
minT = Ts(cellType,1);
maxT = Ts(cellType,2);

spacing = .05;

axisLabelFontSize = 12;
axisTickLabelFontSize = 12;
axisWidth = 1;
izColor = [0 0 0];
glmColor = [.5 .5 .5];


%%% plot filters

f=figure;

subplot('position',[.72+spacing/2 .5+spacing/2 .25-spacing .35-spacing]); hold on;
h1 = gca;
plot([0 length(k)],[0 0],'k--','linewidth',1.5);
plot(k,'linewidth',2);
set(gca,'xtick',0:length(k)/4:length(k),'xticklabel',round(-length(k)*dt:length(k)/4*dt:0))
set(gca,'tickdir','out','linewidth',axisWidth,'fontsize',axisTickLabelFontSize)
text(length(k)/15,max(k)-.05*(max(k)-min(k)),['\mu = ' num2str(round(dc*10)/10)],'fontsize',axisLabelFontSize);
xlim([0 length(k)])
ylim([min(k)-.05*(max(k)-min(k)) max(k)+.05*(max(k)-min(k))])
box off
h1p = get(h1,'position');
set(h1,'position',[h1p(1) h1p(2)*1.05 h1p(3)*.9 h1p(4)*.95])

subplot('position',[.72+spacing/2 .15+spacing/2 .25-spacing .35-spacing]); hold on;
plot([0 length(h)],[0 0],'k--','linewidth',1.5);
plot(h(2:end),'r','linewidth',2);
h2 = gca;
set(gca,'tickdir','out','xtick',0:length(h)/4:length(h),'xticklabel',round(0:length(h)/4*dt:length(h)*dt))
set(gca,'linewidth',axisWidth,'fontsize',axisTickLabelFontSize)
xlabel('time (ms)','fontsize',axisLabelFontSize)
xlim([0 length(h)])
ylim([min(h)-.05*(max(h)-min(h)) max(h)+.05*(max(h)-min(h))])
box off
h2p = get(h2,'position');
set(h2,'position',[h2p(1) h2p(2)*.95 h2p(3)*.9 h2p(4)*.95])

%%% plot simulated data vs. izhikevich data
if minT == 0
    minT = 1;
end
minT = minT/dt;
maxT = maxT/dt;

tIdx = minT:maxT;
t = (tIdx-minT)*dt;

xticks = 0:round(maxT*dt)/4:maxT*dt;

%%% stimulus
subplot('position',[.1 .1+.18*4+spacing/2 .6-spacing .18-spacing]); hold on;
plot(t,I(tIdx),'linewidth',2)
xlim([min(t) max(t)])
ylim([min(I(tIdx))-.05*abs(min(I(tIdx))) max(I(tIdx))+.05*abs(max(I(tIdx)))])
set(gca,'tickdir','out','linewidth',axisWidth,'fontsize',axisTickLabelFontSize,'xtick',xticks,'xticklabel',[])
box off

%%% Izhikevich response
subplot('position',[.1 .1+.18*3+spacing/2 .6-spacing .18-spacing]); hold on;
plot(t,v(tIdx),'color',izColor,'linewidth',1.5)
xlim([min(t) max(t)])
ylim([min(v(tIdx))*1.05 max(v(tIdx))]*1.05)
set(gca,'tickdir','out','linewidth',axisWidth,'fontsize',axisTickLabelFontSize,'xtick',xticks,'xticklabel',[],'ytick',[-80:50:20])
box off

%%% filter outputs ("currents")
Iinj = stimcurr + dc;
subplot('position',[.1 .1+.18*2+spacing/2 .6-spacing .18-spacing]); hold on;
plot(t,Iinj(tIdx),'linewidth',1.5);
plot(t,hcurr(tIdx,1),'r','linewidth',1.5)
xlim([min(t) max(t)])
ylim([min([hcurr(tIdx,1); Iinj(tIdx)])*1.1 max([hcurr(tIdx,1); Iinj(tIdx)])*1.1])
set(gca,'tickdir','out','linewidth',axisWidth,'fontsize',axisTickLabelFontSize,'xtick',xticks,'xticklabel',[])
box off

%%% prob(firing)
subplot('position',[.1 .1+.18+spacing/2 .6-spacing .18-spacing]);
semilogy(t,feval(NL,hcurr(tIdx,1)+Iinj(tIdx)),'color',glmColor,'linewidth',1.5)
xlim([min(t) max(t)])
ylim([.1 10^6])
set(gca,'tickdir','out','linewidth',axisWidth,'fontsize',axisTickLabelFontSize,'ytick',10.^[0 6 12],'xtick',xticks,'xticklabel',[])
box off

%%% GLM spikes
spikeHeight = .7;
subplot('position',[.1 .1+spacing/2 .6-spacing .18-spacing]); hold on;

for i = 1:size(y,2) % for each run of glm simulation
    spt = find(y(tIdx,i));
    for spikeNum = 1:length(spt)
        plot([spt(spikeNum)*dt spt(spikeNum)*dt],[i-.5 i-.5+spikeHeight],'color',glmColor,'linewidth',1.25)
    end
end

spt = find(spikes(tIdx));

for spikeNum = 1:length(spt)
    plot([spt(spikeNum)*dt spt(spikeNum)*dt],[i+.5 i+.5+spikeHeight],'color',izColor,'linewidth',1.25)
    hold on;
end

xlim([0 max(t)-min(t)])
ylim([0 size(y,2)+1+spikeHeight])
set(gca,'tickdir','out','linewidth',axisWidth,'fontsize',axisTickLabelFontSize,'xtick',xticks)
xlabel('time (ms)','fontsize',axisLabelFontSize)


