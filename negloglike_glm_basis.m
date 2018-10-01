function [negloglike,dL,H] = negloglike_glm_basis(prs,NL,xconvki,yconvhi,y,dt,refreshRate,L2pen)

% function [negloglike,dL,H] = negloglike_glm_basis(prs,NL,xconvki,yconvhi,y,dt,refreshRate)
%   prs:        vector of parameters, coefficients for basis functions in order
%                     [kprs (stim filter); hprs (post-spike filter); dc]
%   NL:         function handle for nonlinearity
%   xconvki:    stimulus convolved with each filter basis vector,
%                     upsampled to match response sampling rate
%   yconvhi:    response vector convolved with each filter basis vector
%   y:          response vector (zeros and ones)
%   dt:         time scale of y (in frames/stimulus frame)
%   refreshRate: refresh rate of stimulus (frames/sec)
%   L2pen:      penalty on L2 norm of prs vector
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

negloglike = -y'*g + dt*sum(lambda)/refreshRate + L2pen*prs'*prs;  % negative log likelihood

%% calculate negative gradient
if nargout>1
    
    %%% original code, slower
    %     dL = zeros(size(prs));
    %
    %     % kterms
    %     trm1 = sum(xconvki(logical(y),:),1)';  % sum of above term at t=spike
    %     trm2 = sum(repmat(lambda,1,nkbasis).*xconvki,1)'*dt/refreshRate;
    %     dL(1:nkbasis) = trm2 - trm1; % negative gradient
    %
    %     % hterms
    %     trm1 = sum(yconvhi(logical(y),:),1)';
    %     trm2 = sum(repmat(lambda,1,length(prs)-nkbasis-1).*yconvhi,1)'*dt/refreshRate;
    %     dL(nkbasis+1:end-1) = trm2 - trm1; % negative gradient
    %
    %     % dc
    %     dL(end) = sum(lambda)*dt/refreshRate - sum(y);
    
    dL = zeros(size(prs));
    prsMat = [xconvki yconvhi ones(size(xconvki,1),1)];
    for pr = 1:length(prs)
        dL(pr) = -sum(prsMat(logical(y),pr)) + dt/refreshRate*sum(prsMat(:,pr).*lambda) + L2pen*2*prs(pr);
    end
end

%% calculate negative Hessian
if nargout>2
    
    H = zeros(length(prs));
    
    prsMat = [xconvki yconvhi ones(size(xconvki,1),1)];
    for pr1 = 1:length(prs)
        for pr2 = pr1:length(prs)
            H(pr1,pr2) = dt/refreshRate*sum(prsMat(:,pr1).*prsMat(:,pr2).*lambda) + L2pen*2*(pr1==pr2);
            H(pr2,pr1) = H(pr1,pr2);
        end
    end
    
end
