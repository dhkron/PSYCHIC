function [model, loglikHist] = mixModelFit(data, nmix, type, varargin)
%% Fit a mixture model via MLE/MAP (using EM)
%
% By default we lightly regularize the parameters, so we are doing map
% estimation. To turn this off, set 'prior' and 'mixPrior to 'none'. See
% Inputs below.
%
%% Inputs
%
% data     - data(i, :) is the ith case, i.e. data is of size n-by-d
% nmix     - the number of mixture components to use
% type     - a string, either 'gauss', 'student', or 'discrete'
%            note that 'discrete' means a product of discrete distributions,
%            hence data can (optionally) be vector valued.
% If discrete, we require the observed features to have the same
% number of states/ levels.
%
%% Optional named Inputs
%
% initParams  - a optional struct storing parameters to use for
%               initialization for the first run, (i.e. first
%               randomRestart). The fields depend on the type. For 'gauss',
%               set mu (d-by-nmix), Sigma (d-by-d-by-nmix), mixWeight
%               (1-by-nmix). For 'student', add dof (1-by-nmix) in addition
%               to the gauss' fields. For 'discrete', include T
%              (nObsStates-by-nstates-by-d) as well as mixWeight.
%               If not specified, we initialize 'gauss', and 'student'
%               using kmeans, and 'discrete' by fitting each component on
%               random subsets of the data.
%
% prior      - a struct with type-dependent fields. For 'guass', this is a
%              gaussInvWishart distribution, hence the fields are 'mu', 'Sigma',
%              'dof', 'k'. We don't currently support a prior for
%              'student'. For 'discrete', use 'alpha', the pseudo
%              counts. Alpha must be either a scalar or of size
%              nObsStates-by-nstates-by-d. An alpha value of 1 effectively
%              means no prior since the map estimate adds alpha-1.
%
%              By default, we lightly regularize the parameters, set prior
%              to 'none' to do mle.
%
% mixPrior   - these are pseudoCounts for the mixture distribution. This
%              must be a scalar or of size 1-by-nmix. Set to 'none' to do
%              MLE. A value of 1 effectively means no prior, since the map
%              estimate adds mixPrior -1.
%
% overRelaxFactor - Currently only supported for 'gauss'. If set,
%                   over-relaxed EM is run. See emAlgoAdaptiveOverRelaxed
%                   for details.
%
% See emAlgo for additional EM related inputs
%% Outputs
% 
% A structure - see mixModelCreate for field descriptions
% loglikHist  - a record of the log likelihood at each EM iteration. 
%% 

% This file is from pmtk3.googlecode.com

fprintf('mixModelFit has been deprecated, 2 April 2011\n')

[initParams, prior, mixPrior, overRelaxFactor, EMargs] = ...
    process_options(varargin, ...
    'initParams'        , [], ...
    'prior'             , [], ...
    'mixPrior'          , [], ...
    'overRelaxFactor'   , []);
[n, d]      = size(data);
model.type  = type;
model.nmix  = nmix;
model.d     = d;
mstepOrFn   = [];
model       = setMixPrior(model, mixPrior);
switch lower(type)
    case 'gauss'
        if ~isempty(overRelaxFactor); 
            mstepOrFn = @mstepGaussOr;  
        end
        initFn = @initGauss;
    case 'discrete', 
        initFn = @initDiscrete;
    case 'student'
        initFn = @initStudent;
end
initFn = @(m, X, r)initFn(m, X, r, initParams, prior); 
[model, loglikHist] = emAlgo(model, data, initFn, @estep, @mstep , ...
                            'mstepOR'         , mstepOrFn        , ...
                            'overRelaxFactor' , overRelaxFactor  , ...
                            'verbose', true, ...
                            EMargs{:});
end

%% Initialization
function model = initGauss(model, X, restartNum, initParams, prior)
%% Initialize 
nmix = model.nmix; 
if restartNum == 1
    if ~isempty(initParams)
        mu              = initParams.mu;
        Sigma           = initParams.Sigma;
        model.mixWeight = initParams.mixWeight;
    else
        [mu, Sigma, model.mixWeight] = kmeansInitMixGauss(X, nmix);
    end
else
    mu              = randn(d, nmix);
    regularizer     = 2; 
    Sigma           = stackedRandpd(d, nmix, regularizer); 
    model.mixWeight = normalize(rand(1, nmix) + regularizer); 
end
model.cpd = condGaussCpdCreate(mu, Sigma, 'prior', prior);
end

function model = initDiscrete(model, X, restartNum, initParams, prior)
%% Initialize
nObsStates = max(nunique(X(:)));
if restartNum == 1 && ~isempty(initParams)
    T = initParams.T;
    model.mixWeight = initParams.mixWeight;
else
    % randomly partition data, fit each partition separately, add noise.
    nmix    = model.nmix;
    d       = size(X, 2);
    T       = zeros(nmix, nObsStates, d);
    Xsplit  = randsplit(X, nmix);
    for k=1:nmix
        m = discreteFit(Xsplit{k}, 1, nObsStates);
        T(k, :, :) = m.T;
    end
    T               = normalize(T + 0.2*rand(size(T)), 2); % add noise
    model.mixWeight = normalize(10*rand(1, nmix) + ones(1, nmix));
end
model.cpd = condDiscreteProdCpdCreate(T, 'prior', prior); 
end

function model = initStudent(model, X, restartNum, initParams, prior)
%% Initialize
% Initialize of the back of initGauss, and then just add dof
model = initGauss(model, X, restartNum, initParams, prior); 
if ~isempty(initParams)
    dof = initParam.dof;
else
    dof = 10*ones(1, model.nmix); 
end
initCpd   = model.cpd; 
model.cpd = condStudentCpdCreate(initCpd.mu, initCpd.Sigma, dof, 'prior', prior); 
end

function [ess, loglik] = estep(model, data)
%% Compute the expected sufficient statistics
[weights, ll] = mixModelInferLatent(model, data); 
cpd           = model.cpd;
ess           = cpd.essFn(cpd, data, weights); 
ess.weights   = weights; % useful for plottings
loglik        = sum(ll) + cpd.logPriorFn(cpd) + model.mixPriorFn(model); 
end

function model = mstep(model, ess)
%% Maximize
cpd             = model.cpd;
model.cpd       = cpd.fitFnEss(cpd, ess); 
model.mixWeight = normalize(ess.wsum + model.mixPrior - 1); 
end

%% MstepOR
function [model, valid] = mstepGaussOr(model, modelBO, eta)
%% Over-relaxed Gaussian mstep 
[D, D, nmix] = size(modelBO.cpd.Sigma);
% Since weights are constrained to sum to one, we do update in softmax
% parameterization
mixWeight = model.mixWeight.*(modelBO.mixWeight ./ model.mixWeight).^eta;
mixWeight = normalize(mixWeight);
Sigma     = zeros(D, D, nmix);
mu        = zeros(D, nmix);
valid = true;

muReg    = model.cpd.mu;
muBO     = modelBO.cpd.mu;
SigmaReg = model.cpd.Sigma;
SigmaBO  = modelBO.cpd.Sigma; 

for c = 1:nmix
    % Regular update
    mu(:, c) = muReg(:, c) + eta*(muBO(:, c) - muReg(:, c));
    %Since Sigma is constrained to positive definite matrices, the updation
    %of Sigma is done in the Log-Euclidean space. (ref: "Fast and Simple
    %Calculus on Tensors in the Log-Euclidean Framework", Vincent Arsigny,
    %Pierre Fillard, Xavier Pennec, and Nicholas Ayache)
    try
        matLogSigma    = logm(SigmaReg(:, :, c));
        matLogSigma_BO = logm(SigmaBO(:, :, c));
        matLogSigma    = matLogSigma + eta*(matLogSigma_BO - matLogSigma);
        Sigma(:, :, c) = expm(matLogSigma);
    catch %#ok
        valid  = false; return;
    end
    if ~isposdef(Sigma(:, :, c))
        valid = false; return;
    end
end
model.cpd.mu    = mu;
model.cpd.Sigma = Sigma;
model.mixWeight = mixWeight;
end

