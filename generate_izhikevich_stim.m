function [I, dt] = generate_izhikevich_stim(cellType,T)
%  [I dt] = generate_izhikevich_stim(cellType,T)
%
%  This code generates a stimulus appropriate to the specified type of
%  Izhikevich neuron.  Defaults for each cell type are those used in
%  "Capturing the dynamical repertoire of single neurons with GLMs" (Weber
%  & Pillow 2017).  This code is intended to be used together with
%  "simulate_izhikevich" to produce example responses of Izhikevich
%  neurons.
%
%  Note that the behavior of each Izhikevich neuron type is often not 
%  consistent across different stimulus parameters.  For example, bistable 
%  neurons will only show bistable activity for very particular step times
%  and integration windows (dt).
% 
%  The inputs are:
%   cellType: type of Izhikevich neuron.  
%     Choose from:
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
%   T: max time of stimulus, in ms
%   
%  The outputs are:
%     I: stimulus (current input)
%     dt: time step, in ms  

%% check for unavailable behaviors
if sum([10 17]==cellType)
    disp('This is a subthreshold behavior, which can''t be captured by a GLM.  No example stimulus has been designed for this cell type.')
    I = [];
    dt = [];
    return
end

%% set defaults
if ~exist('T','var') || isempty(T)
    T = 10000;
end


%% parameters
% numbered as in Izhikevich 2004
% a-d values taken from code published by Eugene Izhikevich, found at:
% http://www.izhikevich.org/publications/izhikevich.m
% only 1-9, 11-16, 18-21 make sense for GLM


%       I        dt
pars = [14       0.1 ;    % 1. tonic spiking
        .5       0.1 ;    % 2. phasic spiking
        10       0.1 ;    % 3. tonic bursting
        .6       0.1 ;    % 4. phasic bursting
        10       0.1 ;    % 5. mixed mode
        20       0.1 ;    % 6. spike frequency adaptation
        25       0.1 ;    % 7. Class 1
        .5       0.1 ;    % 8. Class 2
        3.49     0.1 ;    % 9. spike latency
        0        1   ;    % 10. subthreshold oscillations
        .3       0.5 ;    % 11. resonator
        27.4     0.5 ;    % 12. integrator
        -5       0.1 ;    % 13. rebound spike
        -5       0.1 ;    % 14. rebound burst
        2.3      1   ;    % 15. threshold variability
        26.1     0.05;    % 16. bistability
        0        0.1 ;    % 17. depolarizing after-potential
        20       0.1 ;    % 18. accomodation
        70       0.1 ;    % 19. inhibition-induced spiking
        70       0.1 ;    % 20. inhibition-induced bursting
        26.1     0.05 ];  % 21. bistability 2 (Not in original Izhikevich paper)

Ival = pars(cellType,1);
dt = pars(cellType,2);

t = dt:dt:T;

%% generate stimulus
I = zeros(length(t),1);

stepLength = 500; % in units of ms
nStepsUp = floor(T/stepLength/2);

if sum([1:6 10 19 20] == cellType)
    if sum([19 20] == cellType)
        I = 80*ones(length(t),1);
    end
    for i = 1:nStepsUp
        idx = t>(stepLength+stepLength*2*(i-1)) & t<(stepLength*2*(i)+1);
        I(idx) = Ival;
    end
    
elseif sum([7 8] == cellType)
    if cellType == 7
        stepSizes = 15:1:30;
    elseif cellType ==8
        stepSizes = .1:.025:.7;
    end
    for i = 1:length(stepSizes)
        idx = t>(stepLength+stepLength*2*(i-1)) & t<(stepLength*2*(i)+1);
        I(idx) = stepSizes(i);
    end
    
elseif sum([9] == cellType)
    stepLength = 150; 
    nStepsUp = floor(T/stepLength/2);
    
    for i = 1:nStepsUp
        idx = t>(stepLength*1.94+stepLength*2*(i-1)) & t<(stepLength*2*(i)+1);
        I(idx) = Ival;
    end
        
elseif sum([11] == cellType)  % resonator
    stepLength = 150; 
    nStepsUp = floor(T/stepLength/2);
    
    for i = 2:nStepsUp       
        pulseLength = round(5/dt);
        idx = t>(stepLength+stepLength*2*(i-1)) & t<(stepLength+stepLength*2*(i-1)+pulseLength);
        I(idx) = Ival;
     % second pulse
        idx = t>(stepLength+stepLength*2*(i-1)+pulseLength+2*i+pulseLength/2) & t<(stepLength+stepLength*2*(i-1)+2*pulseLength+2*i+pulseLength/2);
        I(idx) = Ival;
    end
elseif sum([12] == cellType)  % integrator
    stepLength = 250; 
    nStepsUp = floor(T/stepLength/2);
    
    for i = 3:nStepsUp
        pulseLength = round(4/dt);
        idx = t>(stepLength+stepLength*2*(i-1)) & t<(stepLength+stepLength*2*(i-1)+pulseLength);
        I(idx) = Ival;
     % second pulse
        idx = t>(stepLength+stepLength*2*(i-1)+pulseLength+6*i+pulseLength/2) & t<(stepLength+stepLength*2*(i-1)+2*pulseLength+6*i+pulseLength/2);
        I(idx) = Ival;
    end
    
    
    
elseif sum([13 14] == cellType)  % rebound spike, rebound burst
    for i = 1:nStepsUp
        idx = t>(stepLength*1.6+stepLength*2*(i-1)) & t<(stepLength*2*(i)+1);
        I(idx) = Ival;
    end
    
elseif sum([15] == cellType)  % threshold variability, steps at random times
    dur = 1/dt; % duration of step in ms
    for i = 1:nStepsUp*2
        idx = stepLength*i-dur:stepLength*i;
        I(idx) = Ival;
        if mod(i,2)
            I(idx-25) = -Ival;
        end
    end
    
elseif sum([16 21] == cellType)
    if cellType == 16
        pulsePolarity = 1;
    elseif cellType == 21
        pulsePolarity = -1;
    end
    
    stepLength = 50;
    nStepsUp = floor(T/stepLength);
    I = I-65;
    pulseDir = 2;
    delay = -3;
    for i = 1:nStepsUp
    
        if mod(i,2)
            idx = t>(stepLength+stepLength*(i-1)) & t<(stepLength+stepLength*(i-1)+pulseDir);
            I(idx) = I(idx)+Ival;
        else
            idx = t>(delay+stepLength+stepLength*(i-1)) & t<(delay+stepLength+stepLength*(i-1)+pulseDir);
            I(idx) = I(idx)+Ival*pulsePolarity;
            
        end
    end
    
elseif sum([18] == cellType)
    baseline = -70;
    I = baseline*ones(size(I));
    for i = 1:nStepsUp
        if mod(i,2)
            idx = t>(stepLength+stepLength*2*(i-1)) & t<(stepLength*2*(i)+1);
            I(idx) = linspace(baseline,baseline+Ival,sum(idx));
        else
            idx = t>(stepLength*1.9+stepLength*2*(i-1)) & t<(stepLength*2*(i)+1);
            I(idx) = linspace(baseline,baseline+Ival,sum(idx));
            
        end
    end
        
end

I = I(1:length(t));
