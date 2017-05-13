function [negloglike,dL] = negloglike_glm_basis_softRect(prs,NL,xconvki,yconvhi,y,dt,refreshRate)

% function negloglike = glm_kh_basis(prs,x,y,dt,kbasis,hbasis)
%   prs:        vector of parameters, coefficients for basis functions in order
%                     [kprs (stim filter); hprs (post-spike filter); dc]
%   NL:         function handle for nonlinearity
%   xconvki:    stimulus convolved with each filter basis vector,
%                     upsampled to match response sampling rate
%   yconvhi:    response vector convolved with each filter basis vector
%   y:          response vector (zeros and ones)
%   dt:         time scale of y (in frames/stimulus frame)
%   refreshRate: refresh rate of stimulus (frames/sec)
%   
% gives same output as negloglike_glm_basis_old, neglogli_GLM, and Loss_GLM_logli
 
%% calculate negative log likelihood
nkbasis = size(xconvki,2); % number of basis functions for k

kprs = prs(1:nkbasis); % k basis functions weighted by given parameters
hprs = prs(nkbasis+1:end-1); %  basis functions weighted by given parameters
dc = prs(end); % dc current (accounts for mean spike rate)

xconvk_dc = xconvki*kprs + dc;

yconvh = yconvhi*hprs; % same as output of spikeconv_mex (line 56, negloglig_GLM)

g = xconvk_dc+yconvh; % g = loglambda for NL = @exp, same as Iinj (line 56, neglogli_GLM)
lambda = feval(NL,g);

% lambda(lambda==Inf) = max(lambda(lambda~=Inf));  % deal with Inf values
loglambda = log(lambda);
% loglambda(loglambda == -Inf) = g(loglambda == -Inf); % deal with values of lambda = 0

negloglike = -y'*loglambda + dt*sum(lambda)/refreshRate;  % negative log likelihood

%% calculate negative gradient
if nargout>1
    
    dL = zeros(size(prs));
    
    % kterms
    trm1all = repmat(exp(g)./lambda./exp(lambda),1,size(xconvki,2)).*xconvki;
    trm1 = trm1all'*y;  % sum of above term at t=spike
    
    trm2all = repmat(exp(g)./exp(lambda),1,size(xconvki,2)).*xconvki;
    trm2 = sum(trm2all)*dt/refreshRate;
    dL(1:nkbasis) = trm2' - trm1; % negative gradient
    
    % hterms
    trm1all = repmat(exp(g)./lambda./exp(lambda),1,size(yconvhi,2)).*yconvhi;
    trm1 = trm1all'*y;  % sum of above term at t=spike
    
    trm2all = repmat(exp(g)./exp(lambda),1,size(yconvhi,2)).*yconvhi;
    trm2 = sum(trm2all)*dt/refreshRate;
    dL(nkbasis+1:end-1) = trm2' - trm1; % negative gradient
    
    % dc
    dL(end) = -y'*(exp(g)./lambda./exp(lambda)) + sum(exp(g)./exp(lambda))*dt/refreshRate; % negative gradient
end
