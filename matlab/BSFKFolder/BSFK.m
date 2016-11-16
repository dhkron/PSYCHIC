function [pp ier] = BSFK(x, y, k, nknots, fixknots, options)
%
% pp = BSFK(x, y)
% [pp ier] = BSFK(x, y, k, nknots, fixknots, options)
%
% Least-squares fitting with Free-Knot B-Spline
%
%     s(x) minimizes sum [s(xi)-y(i)]^2
%          where s(x) is a spline function
%
% INPUTS:
%   Only the first two inputs must be required, the rest is optional
%   - x: abscissa data,
%   - y: data function of x to be fitted by least-squares,
%   - k-1: order of the spline (k=[4] for cubic - defaut; k=2 for linear)
%     the spline (polynomial order on each sub-interval)
%     The spline function has the continuous pth detivative, p=0,1,...,k-2.
%   - nknots: [20] array or scalar. The number of subintervals between two
%     consecutive fixed knots. nknots(:) must be greater or equal to 1.
%     + for array, length(nknots) must be equal to (1 + number of interior
%       fixed-knots).
%     + when nknots is a scalar, it fixes the total number of subintervals.
%       the knots (free and fixed) are distributed more or less
%       equidistance on the whole abscissa interval.
%   - fixknots, position of fixed-knots. Usualy empty array [] (default) is
%     used when user wants all interior knots to be free to move. They must
%     be enclosed by the abscissa data x. 
%   - options: optional structure with following fields:
%     - GPflag: Gauss-newton approximation flag:
%           true: Golub/Pereyra [default],
%           false: Kaufman simplification.
%       Golub/Pereyra formulation is more costly for each iteration.
%       However in some difficult cases, Kaufman formulation might need
%       few more number of iterations to converge, but in most of the
%       configurations Kaufman's is faster.
%     - startingknots: array of starting knots (nknots parameter will be
%       then ignored if startingknots is provided) or either strings
%       'chebyschev' or - default value - 'equidistance'.
%
%   The following parameters are used to control the Gauss-newton stopping
%   criteria:
%
%     - maxiter: maximum number of Gauss-Newton iterations, [100] by
%                default.
%     - sigma: The fit consider to be good when the RMS of the residual
%              is smaller than sigma. By default sigma is estimated
%              automatically from the data.
%       Typically it is recommended to set sigma to the standard deviation
%       of the noise in the signal "y" if it is known.
%       If sigma is not provided, BSFK try to guess the value using the
%       tail of the discrete fourier transform of y and median Chi-square
%       distribution.
%     - epsg: iteration stops when RMS of fit residual changes less than
%                         epsg*sigma. Default is [1e-2].
%     - epsgl: iteration stops when the optimal linsearch solution
%              changes less than epsgl*sigma. Default is [1e-2].
%     - epsk: knots precision; iteration stops when the all the free knots
%             change less epsk*(max(x)-min(x)). Default is [1e-3].
%
%   - display: [0] -> no information will bee displayed on command window.
%                     This is the default value.
%              1 -> fitting progression is traced on command window.
%   - qpengine: engine for quadratic programming. Curently three QP
%     engines are supported
%       1. 'quadprog': Mathworks optimization toolbox is required, or
%       2. 'qpas': by Dr. Adrian Wills, QPC package available at
%                  http://sigpromu.org/quadprog/index.html
%           Bruno strongly recommends using this engine for BSFK.
%       3. 'minq': Prof. Dr. Arnold Neumaier, available from
%                  http://www.mat.univie.ac.at/~neum/software/minq/
%           This engine is NOT supported for shape preserving splines
%           and pointwise constraints
%   - KnotRemoval: 
%       - 'none', no knot-removal step will be carried out; The other
%         values select the strategy used to detect redundant knots;
%       - 'lychemorken' or 'yes' (default),
%       - 'schwetlickschutze' - not working well in this implementation.
%   - knotremoval_factor: scale factor (>1) used to estimate threshold 
%                         for knot-removal. Default value is 1.125.
%     One (1) will be automatically added if its value is smaller than 1.
%     Smaller value -> less knots will be removed, and vice versa.
%   - animation: [0] -> no graphic animation, 1 with animation.
%   - figure: figure handle where the animation is carried out
%             New figure will be created if figure field is not set.
%   - waitbar: boolean [0] or a string of the title of the waitbar
%   - d: derivative-order of the seminorm regularization, d = min(2,k-1)
%        by default.
%   - lambda: regularization parameter, [0] (no regularization) by
%     default. The least-squares is regularized by using the derivative
%     of the spline function normalized on abscissa interval of (0,1):
%     s = argmin J(s) = 1/2 sum (s(xi) - yi)^2 + lambda/2 |R(s)|^2, where
%          R(s) := |d^d(shat)/dx^d|^2,
%               := integral_(0,1) [d^d shat(xi) / dxi^d ]^2 dxi, (eqt1)
%          shat := s [ x/(max(x)-min(x)) ] (normalized spline)
%   - regmethod: 'continuous' (default) or 'discrete'.
%     When regmethod is 'continuous', the exact Gram's matrix defined in
%     (eqt1) is used. This regularization is expensive to compute. When
%     regmethod is set to 'discrete' an equivalent fast discrete version
%     of the regularization is used in place of the continuous (however the
%     number of iterations might increase).
%   - periodic: boolean flag, [FALSE] by default. If periodic is set to
%     TRUE, the k-1 first derivatives (from 0th to k-2th order) of the left
%     and right end are constraints to be equal.
%     This is called "periodic-spline".
%   - 'shape' is an array of substructures that can be used optionally
%     to control the shape of the fitted spline (often known as "shape
%     preserving splines"). The "shape" of the spline is formulated as
%     following inequality constraints on the p^th derivative of the fitted
%     spline function s(x):
%       lo_(p,i) <= d^p(s)/dx^p <= up_(p,i) on the knot interval [ti,ti+1[.
%            where i = 1, 2,..., nknots.
%     The substructure 'shape' must have following fields
%       + 'p': integer, shape-derivative order of the spline s(x) where the
%              constraints will be applied on,
%       + 'lo': lower bound, scalar or array of size nknots,
%       + 'up': upper bound, scalar or array of size nknots.
%     'shape' might contains more than one sub-structures. All the
%     constraints will be taken into account without any attempt of
%     simplification. It is of user's duty to transform the constraints to
%     the simpler but equivalent form before invoking BSFK.
%   - 'pntcon' (stands for "point-wise constraints") is an array of
%     substructures that can be used to enforce the point-wise constraints
%     on the spline function s(x) and/or its derivative. The PNTCON
%     is formulated as following equality constraints on the p^th
%     derivative of the spline function s(x):
%       d^p(s)/dx^p(x_i) = v_i for a set of abscissas {x_1, x_2, ...}
%     The substructure 'pntcon' must have following fields
%       + 'p': integer, pointwise-derivative order of the spline s(x) where
%              the constraints will be applied on,
%       + 'x': abscissa array. Abscissa coordinates of the points.
%            Value smaller than left-bound (e.g., -inf) will be pushed back
%            to the left bound value. Similar change apply for value larger
%            than right bound (e.g., +inf -> right bound).
%       + 'v': constraint values, scalar or array of size nknots where the
%              values/derivatives must meet.
%
% OUTPUTS:
%   - pp: the fit result returned as MATLAB pp-form (piecewise polynomial).
%   - ier: exit code >0 success, <0 faillure.
%
% Various examples are given in TESTBSFK
% Note: The pp structure contains *double* data regardless the class of
%       the inputs.
%
% Reference: algorithm based on Schwetlick & Schutze
% 1. "Least squares approximation by splines with free knots", BIT Num.
%     Math., 35 (1995);
% 2. "Constrained approximation by splines with free knots", BIT Num.
%    Math., 37 (1997).
% 
% See also: spline, ppval, testbsfk
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% Last update: 07-Jun-2010

% History: 17-Nov-2009: Original
%          18-Nov-2009: 'startingknots' parameter 
%                       Simplified discrete regularization
%          19-Nov-2009: remove NaN data before fitting
%                       change TRY/CATCH ME syntax for better compatibility
%                       Estimate automatic of the noise
%          20-Nov-2009: waitbar
%          04-Dec-2009: shape preserving
%          08-Dec-2009: point-wise constraints
%                       Error in Jacobian formula [Schwetlick & Schutze 97]
%                       has been discovered; modify the BSFN accordingly
%          09-Dec-2009: Fix another bug in the Jacobian in constrained case
%          10-Dec-2009: singular constraints will issue a warning (instead
%                       of an error). Refine the Gauss-Newton direction.
%                       Fix few minor bugs.
%          13-Dec-2009: Correct a bug in UpdateConstraints that did not
%                       update the knot position.
%                       Precasting data to double.
%                       Update more frequently the scalling matrix X
%                       for better accuracy when knots get very close.
%                       Reduce the lagrange tolerance to detect active set
%                       of QPC solver to 1e-6*|g| (was 1e-3)
%          18-Mar-2010: fix a small bug of calling minqdef
%          10-Apr-2010: fix a bug with parsing k and nknots
%                       extend to k=1 (piecewise constant fit)
%          31-May-2010: Periodic boundary condition
%          01-Jun-2010: Remove redundant code between BuildJacobian and
%                       simul
%          07-Jun-2010: more robust final conversion to pp form
%                       Symmetrize (eigs 'sa' need true symmetric matrix)
%          08-Jun-2010: Add a workaround when EIGS/ARPACK fails
%          10-Jun-2010: Change in Bernstein.m to avoid NaN for end knots
%          02-July-2011: bug fixed when starting knots are provided
%          06-July-2012: bug fixed, bad checking knots collision (issorted)

% To do: regulatization by GCV or L-curve
%        Boor/Rice initialization
%        knot reinsertion after the first step
%  - check the feasibility of the shape constraints in ConL2Fit (using LP?)

global BSFK_DISPLAY
global WAITBARH
global BSFK_QPENGINE
% persistent ax1

%% Default parameters
if nargin<3 || isempty(k)
    k = 4; % Cubic B-spline
end

if nargin<5 || isempty(fixknots)
    fixknots = [];
end

if nargin<4 || isempty(nknots)
    nsubint = length(fixknots)+1;
    nknots = ceil(20/nsubint)*ones(1,nsubint);
end

if nargin<6 || isempty(options)
    options = struct();
end

if k<1
    error('BSFK: k (=%d) must be larger or equal to 1', k);
end

% columns, we work with sparse, so only double is supported
x = double(x(:));
y = double(y(:));

% Remove NaNs
keep = ~isnan(x+y);
x = x(keep);
y = y(keep);

% sort data in ascending order
[x is] = sort(x);
y = y(is);

% normalize the absissa in [0,1];
xmin = x(1);
xmax = x(end);

% target RMS of fit residual)
m = size(y,1);

%% Default options

% Gauss-newton approximation flag:
% - true: Golub/Pereyra
% - false: Kaufman
GPflag = getoption(options, 'GPflag', true);

% Max number of Gauss-Newton iteration
maxiter = getoption(options, 'maxiter', 100);

% Interation breaking when RMS of residual ge smaller than sigma
sigma = getoption(options, 'sigma', 'auto');
if ischar(sigma) && strcmpi('auto', sigma)
    sigma = Chi2Estimation(y)*0.5; % divided by 2 as safety margin
end

% gradient precision
epsg = getoption(options, 'epsg', 1);
epsgl = getoption(options, 'epsgl', 1e-2);

% knot precision
epsk = getoption(options, 'epsk', 1e-3);

% Level of display, 0 nothing, 1 basic display
BSFK_DISPLAY = getoption(options, 'display', 0);

% Detect the QP engine
qpengine = getoption(options, 'qpengine', GetQPEngine());
% Global variables
BSFK_QPENGINE = qpengine;

% lychemoke, none or 
KnotRemovalMethod = getoption(options, 'KnotRemoval', 'lychemorken');
% Factor to determine the threshold of knot-removal
knotremoval_factor = getoption(options, 'knotremoval_factor', 9/8);
if knotremoval_factor<1
    knotremoval_factor = max(1+knotremoval_factor, 1);
end

% Graphic flag
graphic = getoption(options, 'annimation', ...
          getoption(options, 'animation', false));
Fig = getoption(options, 'figure', []);
if isempty(Fig)
    Fig = {};
else
    Fig = {Fig};
end

% Open graphic figure and basic plot axes
if graphic
    fig = figure(Fig{:});
    %Fig = {fig}; % plot on fig from now on
    clf(fig);
    if bitand(graphic,2^2) || bitand(graphic,2^1)
        ax2 = subplot(2,2,[3 4],'Parent',fig); % bottom
    else
        ax2 = axes('Parent',fig);
    end
end

% waitbar flag
waitbarflag = getoption(options, 'waitbar', 0);

%% Knots generation
% reshape in row, sorted in ascending order, remove both ends
fixknots = setdiff(reshape(fixknots, 1, []), [xmin xmax]);
if any(fixknots<xmin | fixknots>xmax)
    error('BSFK: fixed knots must be in (%g,%g)', xmin, xmax);
end
xscale = (xmax-xmin);
normfun = @(x) (x - xmin) / xscale;
unnormfun = @(x) xmin + xscale*x;
x = normfun(x);
fixknots = normfun(fixknots);
clsx = class(x);

% Generate the initial knot points
% regular spacing in between fixed-knots
startingknots = getoption(options, 'startingknots', 'equidistance');
nknots = max(nknots,1);
t = GuessKnot(nknots, fixknots, clsx, startingknots, normfun);

% Extend knots with k repeated points on left and right brackets
[t knotidx] = extendknots(t, k);

% indice of interior knots (strict)
% Only those knots will be free
p = 1+k:length(t)-k;
% remove indices of fixed-knots
p = setdiff(p, cumsum(nknots(1:end-1))+k);

%% Regularization
lambda = getoption(options, 'lambda', 0);
if lambda<0
    error('BSFK: regularization parameter lambda=%g must be positive', lambda);
end

% Order of the derivative for the regularization term
d = getoption(options, 'd', min(2,k-1));
if lambda && (d>k-1 || d<0)
    error('BSFK: derivative order %d not valid for S%d splines', d, k-1);
end

regmethod = getoption(options, 'regmethod', 'cont');

%% Shape
noshape = [];
shape = getoption(options, 'shape', noshape);
shape = ScaleShape(shape, xscale, size(t,1)-2*k+1);

nopntcon = [];
pntcon = getoption(options, 'pntcon', nopntcon);
pntcon = ScalePntCon(pntcon, xscale, normfun);

% Periodic spline?
periodic = getoption(options, 'periodic', false);

%% Initialization before lauching the loop
niter = 0; % number of Gauss-newton iterations
nouteriter = 0;
validknots = [];
% Array used to keep track of quantities that can trigger iteration stopping
tracearray = nan(maxiter,4);

if BSFK_DISPLAY>=1
    fprintf('BSFK starts\n');
end

if waitbarflag
    timepredict = InitTimePredict(KnotRemovalMethod);
    if ischar(waitbarflag)
        waitbarstr = waitbarflag;
    else
        waitbarstr = 'BSFK progess...';
    end
    WAITBARH = waitbar(0,waitbarstr,'Name','BSFK');
    waitbarflag = true;    
end

stepbystep = true;

% turn off a anoying warning
nrmest_wstate = warning('off', 'MATLAB:normest:notconverge');

% Loop on knots removal
while true
    
    nouteriter = nouteriter + 1;
   
    % Initialize smoothing structure
    % Initialize shape-constrained matrix and rhs (D*alpha >= LU)
    smoothing = InitPenalization(y, t, k, d, lambda, p, regmethod, ...
                                 knotidx, shape, pntcon, periodic);
        
    % no free knot -> nothing else to do
    if isempty(p)
        if ~isempty(validknots)
            smoothing = validknots;
            t = smoothing.t;
            p = smoothing.p;
            ier = 7;
        else
            ier = 8;
        end
        break
    end
    
    % Build the matrix/rhs of Boor and Rice constrained matrix to avoid
    % coalesing of knots: C*t(p) >= h
    % The constraints matrices do *not* depend on free-knots
    [C h] = BuildCMat(p, t);
    
    % Open figure and plot basic stuffs
    if graphic
       
        cla(ax2);
        yfit = fit(x, y, t, k, smoothing);
        plot(ax2, unnormfun(x), y, '.g');
        hold(ax2,'on');
        graph_hyfit = plot(ax2, unnormfun(x), yfit, 'r');
        yl = ylim(ax2);
        plot(ax2, unnormfun(t(1)*[1 1]), yl, 'b-.');
        plot(ax2, unnormfun(t(end)*[1 1]), yl, 'b-.');
        graph_htp = zeros(size(p));
        for j=reshape(p,1,[])
            graph_htp(j) = plot(ax2, unnormfun(t(j)*[1 1]), yl, 'b-.');
        end
        
        if nouteriter==1
            %fixed knots
            pc = setdiff(k:length(t)-k+1,p);
            fknt = t(pc);
        end
        for j=1:length(fknt)
            plot(ax2, unnormfun(fknt(j)*[1 1]), yl, 'c-');
        end
        vcon = find([pntcon.p]==0);
        for j=1:length(vcon)
            s = pntcon(vcon(j));
            plot(ax2, unnormfun(s.x), s.v, 'bo');
        end
        
        drawnow;
        
    end
    
    told = t;
    f = Inf;
    ier = 6;
    
    if waitbarflag
        timepredict.nouteriter = nouteriter;
        timepredict.innertstart = tic();
        timepredict.niterstart = niter;
        timepredict.niter = niter;
    end
    
    % Gauss-Newton iteration
    while niter < maxiter
        
        % Keep track the total number of Gauss-newton iterations
        niter = niter + 1;
        
        if BSFK_DISPLAY>=1
            fprintf('.'); % dots "..."
        end
        
        [r J trash OK] = BuildJacobian(x, t, k, GPflag, smoothing, ...
                                       shape, pntcon, periodic); %#ok
        if ~OK % coalesing knots!!! can't pursue
            t = told;
            ier = -1;
            break
        end
        
        fold = f;
        f = r.'*r;
        tracearray(niter,1) = sqrt(f/m)/sigma;
        if f < m*sigma^2
            ier = 1;
            break
        end
        tracearray(niter,2) = sqrt(abs(f-fold)) / (sigma*epsg);
        if abs(f-fold) < (sigma*epsg)^2
            ier = 2;
            break
        end
        
        % Solve the reduced problem to find descend direction
        s = GaussNewtonStep(t, p, r, J, C, h);
        
        % Line search
        [gammaopt gradline] = linesearch(x, t, p, s, k, smoothing, ...
                                         shape, pntcon, periodic);
        
        tracearray(niter,3) = abs(gradline) / (sigma*epsgl);
        if gammaopt<=0 % lines-search get stuck
            ier = 3;
            break
        elseif abs(gradline)<=sigma*epsgl
            ier = 4;
            break
        end
        
        told = t;
        ds = gammaopt*s;
        dt = expandfreeknot(ds, t, p);
        t = t + dt;
        
        tracearray(niter,4) = max(abs(ds)) / epsk;
        if max(abs(ds)) < epsk
            ier = 5;
            break
        end
        
        % Update waitbar
        if waitbarflag
            timepredict.niter = niter;
            timepredict = TPredict(timepredict, tracearray);
            try, waitbar(timepredict.percent,WAITBARH); end %#ok
        end
        
        if graphic
            
            yfit = fit(x, y, t, k, smoothing);
            set(graph_hyfit, 'YData', yfit);
            for j=reshape(p,1,[])
                set(graph_htp(j), 'XData', unnormfun(t(j)*[1 1]));
            end
            
            drawnow;
            
            if bitand(graphic,2^2)
                gamma = linspace(0,1,51);
                farr = zeros(size(gamma));
                garr = zeros(size(gamma));
                for npnt=1:length(gamma)
                    [farr(npnt) garr(npnt)] = simul(0, gamma(npnt), x, told, ...
                        p, s, k, smoothing, shape, pntcon, periodic);
                end
                
                try, delete(ax1); end %#ok
                if ~bitand(graphic,2^1)
                    ax1 = subplot(2,2,[1 2],'Parent',fig);
                else
                    ax1 = subplot(2,2,1,'Parent',fig);
                end
                [ax1 h1 h2] = plotyy(ax1, gamma, farr, gamma, garr); %#ok
                ylabel(ax1(1), 'f(\gamma)');
                ylabel(ax1(2), '\nablaf(\gamma)');
                
                hold(ax1(1),'on');
                plot(ax1(1),gammaopt*[1 1], ylim(ax1(1)), 'r-.');
                grid(ax1(2),'on');
                if stepbystep
                    reply = input('''CR'' step-by-step; ''c'' continue: ', 's');
                    if ~isempty(reply) 
                        switch lower(reply)
                            case 'c'
                                stepbystep = false;
                            case 's'
                                save('BSFKdebug.mat', '-mat');
                        end
                    end
                end
            end
            

            
        end % graphic animation
        
    end % while-loop Gauss-Newton iteration
    
    if waitbarflag
        timepredict.Elasedtimes(nouteriter) = toc(timepredict.innertstart);
    end
    
    if bitand(graphic,2^1)
        if ~bitand(graphic,2^2)
            ax3 = subplot(2,2,[1 2],'Parent',fig);
        else
            ax3 = subplot(2,2,2,'Parent',fig);
        end        
        semilogy(ax3, tracearray(1:niter,:));
        hold(ax3, 'on');
        grid(ax3, 'on');
        labels = {'sigma','epsg','espgl','espk'};
        state = warning('off','MATLAB:legend:IgnoringExtraEntries');
        legend(ax3,labels{:},'Location','Best');
        warning(state);
    end
    
    if BSFK_DISPLAY>=1
        fprintf('\n'); % a carried return after dots "..."
    end
    
    % check for eventual case of knot-coalesing
    if ~issorted(t)
        t = told;
    end
    
    % No knot-removal
    if strcmpi(KnotRemovalMethod,'none') || (ier<0)
        break
    end
    
    % Keep track the last valid knot sequence before removed the knots
    firstcall = isempty(validknots);
    [t p nremoved validknots knotidx] = ...
        KnotRemoval(x, y, t, k,  smoothing, p, firstcall, ...
                 knotremoval_factor, validknots, shape, pntcon, ...
                 KnotRemovalMethod);
    % Cannot remove anymore knot                        
    if nremoved==0
        % restore the last valid knot state
        if ~isempty(validknots)
            smoothing = validknots;
            t = smoothing.t;
            p = smoothing.p;
        end
        if BSFK_DISPLAY>=1
            fprintf('- nothing can be removed further\n');
            fprintf('Final number of free knots = %d\n', length(p));
        end        
        break; % the while loop
    end

end % while-loop on knot-removal

% Restore the old state of warning
warning(nrmest_wstate);

if waitbarflag
    try, delete(WAITBARH); end %#ok
end

% Do the very last fitting
[yfit alpha r] = fit(x, y, t, k, smoothing);
if BSFK_DISPLAY>=1
    RMS = sqrt(r.'*r / size(r,1));
    fprintf('Final RMS fit residual = %g\n', RMS);
end
        
% The graphic of the final solution
if graphic    
    cla(ax2);
    plot(ax2, unnormfun(x), y, '.g');
    hold(ax2,'on');
    plot(ax2, unnormfun(x), yfit, 'r');
    yl = ylim(ax2);
    plot(ax2, unnormfun(t(1)*[1 1]), yl, 'b-.');
    plot(ax2, unnormfun(t(end)*[1 1]), yl, 'b-.');
    for j=reshape(p,1,[])
        plot(ax2, unnormfun(t(j)*[1 1]), yl, 'b-.');
    end
    for j=1:length(fknt)
        plot(ax2, unnormfun(fknt(j)*[1 1]), yl, 'c-');
    end
    vcon = find([pntcon.p]==0);
    for j=1:length(vcon)
        s = pntcon(vcon(j));
        plot(ax2, unnormfun(s.x), s.v, 'bo');
    end
    
    drawnow;
end

% Convert the result to MATLAB pp-form
pp = Bspline2pp(struct('t', t, ...
                       'k', k, ...
                       'alpha', alpha, ...
                       'unnormfun', unnormfun));

end % BSFK

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
function timepredict = InitTimePredict(KnotRemovalMethod)

if strcmpi(KnotRemovalMethod,'none')
    tratio = 1;
else
    tratio = [5 2 1];
end
timepredict = struct('nouteriter', 0, ...
    'innertstart', tic(), ...
    'niterstart', 0, ...
    'niter', 0, ...
    'Elasedtimes', 0, ...
    'tratio', tratio, ...
    'telapsed', 0, ...
    'tremain', Inf, ...
    'percent', 0);
end
%%
function timepredict = TPredict(timepredict, tracearray)

niter = timepredict.niter - timepredict.niterstart;
dt = toc(timepredict.innertstart);
nouteriter = timepredict.nouteriter;

if  niter>3
    se = Inf;
    for n=1:size(tracearray,2)
        sstart = max(timepredict.niterstart+1,timepredict.niter-5);
        steps = (sstart:timepredict.niter).';
        P = polyfit(steps, log(tracearray(steps,n)), 1);
        if P(1)<0
            sen = -P(2)/P(1);
            if sen>0 && sen<se
                se = sen;
            end
        end
    end
    teinner = (se - timepredict.niterstart) / niter * dt;
    tInnerRemain = teinner - dt;
else
   if timepredict.nouteriter < length(timepredict.tratio)
       n = timepredict.tratio(nouteriter) / timepredict.tratio(1) * 10;
   else
       n = 3;
   end
   teinner = dt;
   tInnerRemain = max((n-niter+1),0)/niter*dt; 
end

teouter = sum(timepredict.Elasedtimes(1:nouteriter-1));
if timepredict.nouteriter < length(timepredict.tratio)
    if nouteriter>1
        r = sum(timepredict.tratio(nouteriter+1:end)) / ...
            sum(timepredict.tratio(1:nouteriter-1));
        tOuterRemain = teouter * r;
    else
        r = sum(timepredict.tratio(nouteriter+1:end)) / ...
            timepredict.tratio(1);
        tOuterRemain = (teinner + tInnerRemain) * r;
    end
else
    tOuterRemain = 0;
end

timepredict.telapsed = teouter + teinner;
timepredict.tremain = tInnerRemain + tOuterRemain;
timepredict.percent = timepredict.telapsed / ...
                     (timepredict.telapsed + timepredict.tremain);
                 
end

%%
function sigma = Chi2Estimation(y)
% estimate the standard deviation of noise from the data

n = size(y,1);
yfft=fft(y(:)) / sqrt(n);

nremoved = max([1 ceil((n-1)/16)]);
nremoved = min(nremoved,floor((n-1)/2));
yfft([1+(0:nremoved-1) end-(0:nremoved-1)]) = [];
m = median([real(yfft); imag(yfft)].^2);
sigma = sqrt(m/(1-2/3+4/27-8/729));

end

%%
function qpengine = GetQPEngine()
% function qpengine = GetQPEngine()
% Detect the QP engine installed on the computer

if exist('quadprog','file') >= 2
    qpengine = 1;
    return;
end
if exist('qpas','file') == 3 % mex
    qpengine = 2;
    return;
end
qpengine = 3;

end % GetQPEngine

%%
function t = GuessKnot(nknots, fixknots, cls, startingknots, normfun)
% t = GuessKnot(nknots, fixknots, cls, startingknots, normfun)
% Generate the initial knot points
% regular spacing in between fixed-knots

% Number of subintervals (subintervals is delimited by 2 fixed-knots)
nsubint = length(fixknots)+1;
if length(nknots)==1
    fixknots = reshape(fixknots, 1, []);
    a = diff([0 fixknots 1]);
    n = nknots*(a / sum(a));
    [ns is] = sort(n);
    ns(end) = ns(end) + (nknots-sum(ns));
    for i=1:length(n)+1
        ni = [ceil(ns(1:i-1)) floor(ns(i:end))];
        if sum(ni)==nknots            
            break
        end
    end
    nknots(is) = max(ni,1); % make sure there is there is at least
                            % 1 interval between two adjadcent fixed-knots
elseif length(nknots) ~= nsubint
    error('BSFK: incompatible length(nknots) and length(fixknots)+1 (%d,%d)', ...
          length(nknots), nsubint);
end

if nargin<4 || isempty(startingknots)
    startingknots = 'equidistance';
else
    if isnumeric(startingknots)
        startingknots = reshape(startingknots, 1, []); 
        startingknots = normfun(startingknots);
        t = union(startingknots, fixknots);
        % Make sur the knots cover the data
        if t(1)>0
            t = [0 t];
        end
        if t(end)<1+eps
            t(end+1) = 1+eps;
        end
        t = t(:); % Bug fixed, reported by Weathers, Douglas Ervin June 28, 2011
        return
    end
end

istart = [0 cumsum(nknots)]+1;
t = zeros(1,istart(end), cls);
for i=1:length(nknots)
    if i==1
        left = 0;
    else
        left = fixknots(i-1);
    end
    if i==nsubint
        % eps to avoid end effect of the right brackets
        right = max([1+eps left]);
    else
        right = fixknots(i);
    end
    switch startingknots
        case 'chebyschev'
            % Chebyschev knots
            ti = linspace(0, pi, nknots(i)+1);
            ti = (1-cos(ti))/2;
            ti = left + ti*(right-left);
        otherwise
            if ischar(startingknots)
                % equidistance
                ti = linspace(left, right, nknots(i)+1);
            end
    end
    subs = istart(i):istart(i+1)-1;
    t(subs) = ti(1:end-1);
end
t(end) = right;
t = t(:);

end % GuessKnot

%%
function [t knotidx] = extendknots(t, k)
% function [t knotidx] = extendknots(t, k)
%
% Extend the knots by repeating with k times in both ends
%
% knotidx: keeps track of the original index of the knots (before
% removing step, but with knots at two ends repeat k times)

lt = length(t);
idx = [1+zeros(1,k) 2:lt-1 lt+zeros(1,k)];
t = t(idx);
knotidx = 1:length(t);

end % extendknots

%%
function t = expandfreeknot(freeknot, t, p)
% function t = expandfreeknot(freeknot, t, p)
% expand the free-knots to entire set of knots by padding zeros
% for fixed knots
% NOTE: values of T not used, only the size

t = zeros(size(t));
t(p) = freeknot;

end

%%
function smoothing = InitSmoothing(y, t, k, d, lambda, p, ...
                                   regmethod, knotidx, periodic)
% smoothing = InitSmoothing(y, t, k, d, lambda, p, regmethod, knotidx, ...
%                           periodic)
% Initialization the smoothing structure with some quantities we need
% later on.
% knotidx: keeps track of the original index of the knots (before
% removing step, but with knots at two ends repeat k times)

% Number of free nodes
l = length(p);
% total number of knots
nk = length(t);
% Basis of free-knots
E = accumarray([p; 1:l].', 1, [nk l]);

active = lambda>0;
if active
    sqrtlbd = sqrt(lambda);
    q = length(t)-(k+d); % dimension of the regularization mapping
    r = k-d;
    if ischar(regmethod) % convert to numerical code
        switch lower(regmethod)
            case {'cont' 'continuous' 'c' 'gram'}
                regmethod = 1;
            case {'disc' 'discrete' 'd' 'dgd'}
                regmethod = 2;
            otherwise
                error('BSFK: unknown regularization method %s', regmethod);
        end
    end
    knotsubsidx = 1+d:length(t)-d;
    smoothing = struct('active', true, ...
        'regmethod', regmethod, ...
        't', t, ...
        'knotidx', knotidx, ... % original index of the knots 
        'k', k, ...
        'p', p, ...
        'lambda', lambda, ...
        'sqrtlbd', sqrtlbd, ...
        'd', d, ... % order of the derivative to be regularized
        'subsidx', knotsubsidx, ...
        'r', r, ... % polynomial-order of the d-derivative
        'm', size(y,1), ... % dimension of fitting data
        'q', q, ... % dimension of the space of the d-derivative of the spline space
        'ydata', y, ...
        'ytarget', cat(1, y, zeros(q,1)), ...
        'E', E, ...
        'GramE', E(knotsubsidx,:), ... % retrict the direction to relevant knot indexes of Gram matrix
        'D', [], 'LU', [], 'X', [], ... % shape-constraints, update later
        'Deq', [], 'veq', [], 'Xeq', [], ... % point-wise constraints, update later
        'periodic', periodic); % Periodic spline?
else
    smoothing = struct('active', false, ...
        'regmethod', 'none', ...
        't', t, ...
        'knotidx', knotidx, ... % original index of the knots 
        'k', k, ...
        'p', p, ...
        'lambda', 0, ...
        'sqrtlbd', 0, ...
        'd', 0, ...
        'subsidx', 1:length(t), ...
        'r', k, ...
        'm', size(y,1), ...
        'q', 0, ...
        'ydata', y, ...
        'ytarget', y, ...
        'E', E, ...
        'GramE', zeros(length(t),0,class(t)), ...
        'D', [], 'LU', [], 'X', [], ...
        'Deq', [], 'veq', [], 'Xeq', [], ... 
        'periodic', periodic);
end

end % InitSmoothing

%%
function [alpha actset uineq ueq] = ConL2Fit(y, B, D, LU, X, Deq, veq, Xeq, ...
                                     qpengine)
% [alpha actset u] = ConL2Fit(y, B, D, LU, X, Deq, veq, Xeq)
% Call ConL2Fit(..., qpengine) to specify the QP engine
%
% Solve the fixed-knot least-squares fitting problem under
% shape-constraints by using quadratic programming
%
% Minimize J(alpha) := 1/2 | B*alpha - y |^2
%          such that X*(D*alpha -LU) >= 0
%
% alpha: is the optimal solution
% actset: is the set of active constraints, i.e.,
%         X*(D*alpha - LU)_i >= 0, for all i in ACTSET
% u: Lagrange multiplier
% 

global BSFK_DISPLAY
global BSFK_QPENGINE

% By default, use the global QP engine set by the main function
if nargin<9 || isempty(qpengine)
    qpengine = BSFK_QPENGINE;
end

% Hessian and residual gradient vector for QP
H = B.' * B;
g = -(y.' * B).';

cls = class(g);

% Normalized constraints
D = X*D;
LU = X*LU;
Deq = Xeq*Deq;
veq = Xeq*veq;



% Call the QP routine to solve constrained quadratic problems
% Minimize J(alpha) := 1/2 | B*alpha - y |^2
%          such that D*alpha >= LU
switch lower(qpengine)
    
    case {1, 'quadprog'}, % mlint does not like mixed type
        
%         error(['BSFK: quadprog solver is supported for shape-preserving splines\n' ...
%                'QPC engine is recommended']);
        % Blind coded by the author, who doesn't own optimization toolbox
        
        % Tolerance on lagrange parameters to detect active set.
        lagrangetol = 1e-3*norm(g);
        
        qpoptions = optimset('Display','off',...
                             'LargeScale','off');
        s0 = zeros(size(g),cls);
        Aeq = Deq; beq = veq;
        lb = []; ub = [];
        % medium-scale algorithm, because rge presense of linear inequalities
        [alpha trash ier qpout lbd] = quadprog(H, g, -D, -LU, ...
                                  Aeq, beq, lb, ub, s0, qpoptions); %#ok
        if ier<0 && BSFK_DISPLAY>=1
            warning('BSFK:QPNonConvergence', ...
                    'ConL2Fit: not convergence encountered by MINQDEF');
        end
               
        % dual variable
        uineq = lbd.ineqlin(:);
        ueq = lbd.eqlin(:);
        % Check the sign of the lagrange multiplier
        s = (ueq.'*Deq)*(H*alpha+g);
        ueq = sign(s)*ueq;
        
        actset = find(uineq>lagrangetol);
        uineq = X*uineq; % X.'*u; but X is symmetric
        ueq = Xeq*ueq; % idem
        
    case {2, 'qpas' 'qpc'}, % mlint does not like mixed type
        
        % Tolerance on lagrange parameters to detect active set.
        lagrangetol = 1e-6*norm(g);
        
        A = [D; Deq];
        b = [LU; veq];
        eq = [false(size(LU)); true(size(veq))];
        qpoptions = struct('engine','qpas',...
                           'prt',false);
        [alpha dual ier lbd] = qpmin(g, H, A, b, eq, qpoptions); %#ok
        if ier && BSFK_DISPLAY>=1
            warning('BSFK:QPNonConvergence', ...
                    'ConL2Fit: possibly unfeasible shape constraints');  
        end
        
        % dual variable
        uineq = lbd.inequality(:);
        ueq = lbd.equality(:);
        % Check the sign of the lagrange multiplier
        s = (ueq.'*Deq)*(H*alpha+g);
        ueq = sign(s)*ueq;
        
        actset = find(uineq>lagrangetol);
        uineq = X*uineq; % X.'*u; but X is symmetric
        ueq = Xeq*ueq; % idem
        
    otherwise % MINQDEF engine, written in Matlab
        error('BSFK:ShapeQPengine',['BSFK: MINQDEF solver is supported for shape-preserving splines\n' ...
                 'QPC engine installation is recommended\n' ...
                 'from http://sigpromu.org/quadprog/index.html']);
%         A = [D; Deq];
%         b = [LU; veq];
%         eq = [false(size(LU)); true(size(veq))];
%         qpoptions = struct();
%         prt = false;
%         [alpha dual ier]=minqdef(g, H, A, b, eq, prt, qpoptions); %#ok
%         
% %         if ier && BSFK_DISPLAY>=1
% %             warning('BSFK:QPNonConvergence', ...
% %                 'ConL2Fit: not convergence encountered by MINQDEF');
% %         end
% 
%         dual = dual(:);
%         % dual variable
%         uineq = dual(1:size(D,1));
%         ueq = -dual(size(D,1)+1:end);
%         
%         actset = find(uineq>lagrangetol);
%         uineq = X*uineq; % X.'*u; but X is symmetric
%         ueq = Xeq*ueq;
        
end % switch

end % ConL2Fit

%%
function [B out2 alpha yfit r dB] = ModelMat(x, t, k, smoothing, E)
% [B Bi alpha yfit r] = ModelMat(x, t, k, smoothing) OR
% [B Bi alpha yfit r dB] = ModelMat(x, t, k, smoothing) OR
% [B Bi alpha yfit r dB] = ModelMat(x, t, k, smoothing, E)
%
% Compute the model matrix (B),
% its pseudo-inverse (Bi),
% the inversion coefficients (alpha),
% the data fitting (yfit), including regularization target,
% the residual vector (r),
% optional the knot-derivative (dB) in the direction E. If E is not
% provided, all the free-knots directions will be calculated.
% 
% When shape-constraint is active, the active set is returned in the
% ans lagrange second output (and not the pseudo-inverse matrix Bi,
% irrelevant for the next):
%      [B s alpha yfit r dB] = ModelMat(...)
% s is a structure with s.active contains active set
%                       s.lagrange is lagrange parameters

% Do we need to compute the derivative?
isdBNeeded = nargout >= 6;

if isdBNeeded
    if nargin<5 % E is not provided
        E = smoothing.E;
    end
    [B dB] = BernKnotDeriv(x, t, [], k, E);  % 50 ms
else
    B = Bernstein(x, t, [], k);
end

if isempty(B) % in case knots bump into each other (floating point
              % innacuracy + violation of Schoenberg-Whitney condition)
    % Return all remaining outputs as empty
    [out2 alpha yfit r dB] = deal([]);
    return
end

if smoothing.active
    
    % (n+k-d) knots for (n-d) derivative basis
    tr = t(smoothing.subsidx);
    
    if isdBNeeded
        if nargin<5
            % E is already assigned few lines above
            GramDknot = {smoothing.GramE};
        else
            % Restrict the direction to the relevant knots 'tr' for Gram
            GramDknot = {E(smoothing.subsidx,:)};
        end
        DDknot = {E};
    else
        GramDknot = {};
        DDknot = {};
    end

    % d-derivative matrix
    [Dr td kd dDr] = DerivBKnotDeriv(t, k, smoothing.d, DDknot{:}); %#ok
    
    % Compute the L^2 scalar product Gram matrix and the knot-derivatives
    switch smoothing.regmethod
        case 1,
            [R dR] = Gram(tr, smoothing.r, GramDknot{:}); % Gram matrix   
        case 2,
            [R dR] = DGD(tr, smoothing.r, GramDknot{:}); % Discrete matrix   
    end
    
    if isdBNeeded % set the knot derivative-matrix dB
        
        p = size(E,2); % number of derivatives
        % Cholesky factorization and derivative
        [U dU] = derivchol(R, dR);
        
        % Compute the derivative of the product U*Dr
        if p==1 % dDr is not returned in cell by DerivBKnotDeriv
            dUD = dU*Dr + U*dDr;
        else
            dUD = zeros(size(U,1),size(Dr,2),p,class(t));
            % Loop over the directions
            for i=1:p
                dUD(:,:,i) = dU(:,:,i)*Dr + U*dDr{i};
            end
        end
        % Concatenate derivatives together
        dB = cat(1, dB, smoothing.sqrtlbd*dUD);
    else
        % Cholesky factorization
        U = derivchol(R);
    end
    % Concatenate model and regularization
    B = cat(1, B, smoothing.sqrtlbd*U*Dr);
end

y = smoothing.ytarget;
D = smoothing.D;
Deq = smoothing.Deq;
if isempty(D) && isempty(Deq)
    Bi = pseudoinverse(B);
    % Uncontrained regression with fixed knots
    alpha = Bi*y;
    out2 = Bi;
else
    % Shape-constrained regression
    [alpha actset uineq ueq] = ConL2Fit(y, B, D, smoothing.LU, smoothing.X, ...
                                Deq, smoothing.veq, smoothing.Xeq);
    out2 = struct('actset', actset, ...
                  'lagrange_ineq', uineq, ...
                  'lagrange_eq', ueq);
end

% yfit is the fit data, in other word the projection on span(B) of ytarget
yfit = B*alpha;

% residual vector
r = smoothing.ytarget - yfit;

end % ModelMat
               
%%
function [yfit alpha r] = fit(x, y, t, k, smoothing)
% function [yfit alpha r] = fit(x, y, t, k, smoothing)
% Return
%   yfit: current fit to y with knots t
%   alpha: the B-spline coefficients
%   r: residual vector (model alone)

% Update the constraints matrices corresponding to new knot positions
smoothing = UpdateConstraints(smoothing, t);

[B Bi alpha yfit] = ModelMat(x, t, k, smoothing); %#ok
if isempty(B) % knots bump to each-other
    yfit = nan(size(y));
else
    % Remove the regularization trailing
    yfit = yfit(1:size(y,1));
end

r = y-yfit;

end % fit

%%
function [t p nremoved validknots knotidx] = KnotRemoval(x, y, t, k, smoothing, p, ...
                             firstcall, knotremoval_factor, validknots, ...
                             shape, pntcon, strategy)
% [t p rmidx validknots knotidx] = KnotRemoval(x, y, t, k, smoothing, p, ...
%                             firstcall, knotremoval_factor, validknots, ...
%                             shape)
%
% Try the remove the redundant knots.
% INPUTS:
% - (x,y) are data; k is spline order
% - t: current knot-vector
% - p: current vector of index of free-knots
% - smoothing: current smoothing vector
% - firstcall: boolean, must be set to TRUE when KnotRemoval is invoked the
%              first time so that KnotRemoval computes and preserves the
%              residual threshold for later use
% - knotremoval_factor: a scale factor (>1) used to estimate threshold
% - validknots: a current valid knot state (i.e., redisual smaller than
%               the threshold). This is the same structure as SMOOTHING
%               when FIRSTCALL input is true, validknots is not used.
%               KnotRemoval will update the VALIDKNOTS structure and
%               return to the caller
% -shape: shape-constrained structure
% -strategy: 'lychemorken' or 'schwetlickschutze' (method selection)
% OUTPUTS:
% - t: new knot-vector with less knots
% - p: corresponding vector of index of free-knots
% - nremoved: number of knot that has been removed
% - validknots: see descrption in INPUTS
% - knotidx: keeps track of the original index of the knots (before
%            removal, but with knots at two ends repeat k times)

global BSFK_DISPLAY
persistent THRESHOLD

% Lyche & Morken knot-removal strategy
if nargin<12 || isempty(strategy)
    strategy = 'lychemorken';
end

% Current knot index, it will be modified when a knot is removed
knotidx = smoothing.knotidx;

if BSFK_DISPLAY>=1
        fprintf('knot-removing ');
end
    
% Number of knot removed
nremoved = 0;
% Loop trying to remove sequentially knots one-by-one
while true
    
    % Compute the residual to data
    [yfit alpha r] = fit(x, y, t, k, smoothing); %#ok
    rmsresidu = sqrt((r.'*r)/size(r,1));
    
    if firstcall
        % Save the threshold in persistent variable for later use
        THRESHOLD = knotremoval_factor*rmsresidu;
        firstcall = false;
    end
 
    % Oops, residual increases too much when this knot would be removed
    if rmsresidu >= THRESHOLD
        if BSFK_DISPLAY>=1 && nremoved>0
            fprintf('\n%d knots removed\n', nremoved);
        end
        return
    else
        validknots = smoothing;
        validknots.t = t;
    end
    
    if isempty(p)
        % nothing left
        return
    end
    
    if BSFK_DISPLAY>=1 && nremoved>0
        fprintf('-');
    end
    
    switch lower(strategy)
        
        case {'lychemorken', 'yes'} % preferable strategy
            % Find the knot so that the fit residual is smallest
            % after it has been removed
            minrms = Inf;
            % Loop on all knots
            for ir=1:length(p)
                [tr pr sr] = remove(ir); %#ok
                [yfit alpha r] = fit(x, y, tr, k, sr); %#ok
                rmsir = sqrt((r.'*r)/size(r,1));
                if rmsir<minrms % keep track the best
                    minrms = rmsir;
                    iremoved = ir;
                end
            end % for-loop
            
        case 'schwetlickschutze' % Schwetlick & Schutze, does not work very 
            % efficiently, may be BL miss some trick
            % Compute the jump of the first discontinous derivative (k-1)
            tp = t(p);
            Brleft = Bspline(tp, t, k, alpha, k-1, true);
            Brright = Bspline(tp, t, k, alpha, k-1, false);
            JumpBr = abs(Brright-Brleft);
            [trash iremoved] = min(JumpBr); %#ok
            
        otherwise
            
            error('BSFK: unknown knot-removal strategy %s', strategy);
    end

    % Remove the selected knot for good
    [t p smoothing knotidx] = remove(iremoved);
    nremoved = nremoved + 1;
    
end % while loop

    % Nested function return a temporary state when the knot #iremoved
    % is removed
    function [tr pr sr knotidx] = remove(iremoved)
        [tr pr] = deal(t, p); % copy t and p
        tr(p(iremoved)) = [];
        knotidx = smoothing.knotidx;
        knotidx(p(iremoved)) = [];
        pr = [p(1:iremoved-1) pr(iremoved+1:end)-1];
        periodic = smoothing.periodic;
        sr = InitPenalization(y, tr, k, smoothing.d, smoothing.lambda, pr, ...
                              smoothing.regmethod, knotidx, ...
                              shape, pntcon, periodic);
    end

end % KnotRemoval

%%
function shape = ScaleShape(shape, xscale, nknots)
% shape = ScaleShape(shape, xscale)
% Scale the shape lower/upper bounds to accomodate the fact that we change
% the scale of the abscissa x.

if isempty(shape)
    % a structure for no constraint
    shape = struct('p',{}, 'lo', {}, 'up', {});
else
    l = length(shape);
    P = [shape.p];
    if numel(unique(P)) < numel(P)
        warning('BSFK:ShapeDuplicateOrder', 'BSFK: SHAPE has duplicated order');
    end
    for i=1:l
        s = shape(i);
        p = s.p;
        lo = s.lo;
        up = s.up;
        if numel(lo)~=1 && numel(lo)~=nknots
            error('BSFK: LO shape %d is not valid (%d elements expected)', ...
                   i, nknots);
        end
        if numel(up)~=1 && numel(up)~=nknots
            error('BSFK: UP shape %d is not valid (%d elements expected)', ...
                   i, nknots);
        end
        a = xscale^p;
        shape(i).lo = a*lo;
        shape(i).up = a*up;
    end % for-loop
end

end % ScaleShape

%% 
function pntcon = ScalePntCon(pntcon, xscale, normfun)
% pntcon = ScalePntCon(pntcon, xscale)
% Scale the pntcon abscissa/derivative values to accomodate the fact that
% we change the scale of the abscissa x.

if isempty(pntcon)
    % a structure for no constraint
    pntcon = struct('p',{}, 'x', {}, 'v', {}, 'nc', {});
else
    l = length(pntcon);
    P = [pntcon.p];
    if numel(unique(P)) < numel(P)
        warning('BSFK:PntConDuplicateOrder', 'BSFK: PNTCON has duplicated order');
    end
    for i=1:l
        s = pntcon(i);
        p = s.p;
        % linearly map to [0,1]
        x = normfun(s.x(:));
        nc = size(x,1);
        if numel(s.v)~=1 && numel(s.v)~=nc
            error('BSFK: V pntcon %d is not valid (%d elements expected)', ...
                   i, nc);
        end
        a = xscale^p;
        v = zeros(size(x),class(s.v));
        v(:) = s.v;
        % truncated to left/right bounded
        x = max(min(x,1),0);
        % use the last value if there is conflict, x is sorted
        [x is] = unique(x);
        pntcon(i).x = x;
        pntcon(i).v = a*v(is); % scaled derivative
        pntcon(i).nc = nc;
    end % for-loop
end

end % ScalePntCon

%% 
function [D out2 X] = BuildDineqMat(t, knotidx, k, shape, alpha, lagrange, E)
% [D LU X] = BuildDineqMat(t, knotidx, k, shape);
%
% Initialize shape inequality constraint matrix "D" and rhs "LU" such that
%   the constraints: D*alpha >= LU, (alpha beeing the vector of coeffs)
% are sufficient conditions to ensure the desired shape of the spline
%
% Note: D is sparse, LU is full
% shape will contains modified shape with initialization
%
% [D LU X] = BuildDineqMat(...) return a scaling diagonal matrix X such that
% X*D is normalized. The shape-constraints can be equivalently handled
%        (X*D) * alpha >= X*LU
% for better numerical stability
%
% Call with single output to get the update solely of the matrix
%    D = BuildDineqMat(t, knotidx, k, shape);
%
% If alpha, lagrange, and E inputs are provided, the knot-derivative
% is returned and its adjoint
%   [Dt Dtt] = BuildDineqMat(t, knotidx, k, shape, alpha, lagrange, E);
%

global BSFK_DISPLAY

% current number of basis splines
n = size(t,1) - k;

% number of shape constraints
l = length(shape);
P = [shape.p];
if any(P>=k)
    ibad = find(P>=k,1,'first');
    error(['BSFK: shape derivative order p(=%d) ' ...
           'must be smaller than k(=%d)'], P(ibad), k);
end

%% In the following for loop, we will build the lower and upper bounds for
% the coefficients of the derivatives in the B-spline spaces. This is a
% *sufficient* conditions (strong) to enforce the desired shape.
% The result in this step does not depend on the knot locations; it remains
% identical for the same configuration of knots when betwen two occurences
% of knot-removal procedures.

cls = class(t);

if nargout>=2 && nargin<=4 % rhs is requested
    % the number of complete knots (when nothing has been removed)
    norg = knotidx(end);
    nknots = norg+1-2*k; % number of subintervals when BSFK started
    
    % Because of knot-removal precedure, a new subintervals might contain
    % many original subintervals. We run HISTC to which contain which.
    % We start with k, because there is (k-1) duplicated knots on the
    % left side, and there are nknots original interval
    [trash subint] = histc(k:k+nknots-1,knotidx); %#ok
    % put in column
    subint = subint(:);
    lo = zeros(nknots,1,cls);
    up = zeros(nknots,1,cls);
    % Loop on all constraint derivative
    for i=1:l
        p = shape(i).p;
        % number of intervals where the pth-derivative B-basis functions
        % lay on
        w = k-p;
        % first and last intervals in which the the pth-derivative of the
        % spline has its support
        firstint = 1+p;
        lastint = n+k-p-1;
        
        % Take the max of lower bounds where sub-intervals are now grouped
        % together due to knot-removal
        lo(:) = shape(i).lo(:);
        logrouped = accumarray(subint, lo, [lastint 1], @max, -Inf(cls));
        % runing max on logrouped, there will be (n-p) LO
        LO = minmaxfilt(logrouped(firstint:lastint), w, 'max', 'valid');
        % This is the lower bound of the pth derivative coefficients
        shape(i).LO = LO;
        
        % Take the min of upper bounds where sub-intervals are now grouped
        % together due to knot-removal
        up(:) = shape(i).up(:);
        upgrouped = accumarray(subint, up, [lastint 1], @min, +Inf(cls));
        % runing min on upgrouped, there will be (n-p) UP
        UP = minmaxfilt(upgrouped(firstint:lastint), w, 'min', 'valid');
        % This is the upper bound of pth derivative coefficients
        shape(i).UP = UP;
    end % for-loop
    
    %% Build the the rhs LU, full matrix
    % The result in this step does not depend on the knot locations;
    
    % the lower bound will be at put in the first-half of the matrices
    % upper bound in the second-half
    if isempty(shape)
        LU = zeros(0,1);
    else
        LO = cat(1,shape.LO);
        UP = cat(1,shape.UP);
        LU = [+LO; ...
              -UP]; % (2*m) x 1
    end
    
    out2 = LU;
end % if nargout>=2 && nargin<=4

%% Build the linear inequality matrix D
% The result in this step does depend on the knot locations

% half number of constraints
m = sum(n-P);
    
if nargin<=4
    % number of non-zeros elements of D
    % for each structure, there is (n-p) rows
    % each row has (p+1) non zeros elements
    nzD = sum((n-P).*(P+1));
        
    % Loop on all constraint derivative
    % building the indexes and value of the linear constraint matrix
    ilast = 0;
    row = zeros(nzD,1,'double');
    col = zeros(nzD,1,'double');
    val = zeros(nzD,1,'double'); % sparse does not support single 
    xn = zeros(m,1,'double');
    rowstart = 0;
    singtol = eps(1);
    for i=1:l
        % Dp is sparse of the size (n-p) x n
        Dp = DerivB(t, k, shape(i).p); % Derivative matrix of order p^th
        % it should find (n-p)*(p+1) elements
        [r c dp] = find(Dp);
        
        nz = length(dp); % == (n-p)*(p+1)
        mi = size(Dp,1); % == (n-p)
        
        % lower constraints
        r = rowstart+r;
        idx = ilast+(1:nz);
        row(idx) = r;
        col(idx) = c;
        val(idx) = dp;
        dn = normest(Dp);
        if BSFK_DISPLAY && (dn<singtol)
            warning('BSFK:SingularCon', ...
                    'pointwise-constrained matrix close to singular');
        end
        xn(rowstart+(1:mi)) = 1/dn;
        % update the counters
        rowstart = rowstart+mi;
        ilast = ilast + nz;
    end % for-loop for D
    
    if any(isinf(xn))
        error('BSFK: singular pointwise-constrained matrix');
    end
    
    % Build the matrix for lower constraints
    if ilast<nzD
        idx = (1:ilast);
        D = sparse(row(idx),col(idx),val(idx),m,n);
    else % no need for subindexing
        D = sparse(row,col,val,m,n);
    end
    % (2*m) x n
    % the upper constraints are minus the lower constraints
    D = [D; -D];
    
    % Left diagonal scaling matrix (m x m)
    if nargout>=3
        X = spdiags([xn; xn], 0, 2*m, 2*m);
    end
else % Build the derivative of D with respect to knots
    
    % Loop on all constraint derivative
    % building the indexes and value of the linear constraint matrix
    ne = size(E,2);
    rowstart = 0;
    DL = zeros([n 1 ne], 'double');
    DR = zeros([m ne], 'double');
    
    % lagrange(:,1) lower bound lagrange vector
    % lagrange(:,2) upper bound lagrange vector
    lagrange = reshape(lagrange, [], 2);
    ul = lagrange(:,1);
    uu = -lagrange(:,2);
    for i=1:l
        p = shape(i).p;
        % DDp is full of the size (n-p) x size(E,2)
        mi = n-p;
        idx = rowstart+(1:mi);
        uil = ul(idx);
        uiu = uu(idx);
        [alphad td kd Dr trash Dl] = ...
            DerivBKnotDeriv(t, k, p, E, alpha, [uil uiu]); %#ok
        DR(idx,:) = Dr;
        % accumulate the gradient of the constraints
        DL = DL + sum(Dl,2);
        rowstart = rowstart+mi;
    end % for-loop for D
    
    % (2*m) x size(E,2)
    % return the knot-derivative of the matrix as the first output
    D = [DR; -DR];
    
    % n x size(E,2)
    % return the dual knot-derivative (Dtt) of the matrix as the second
    % output
    out2 = reshape(DL, [n ne]);
    
end % if nargin<=4

end % BuildDineqMat

%%
function [Deq out2 Xeq] = BuildDeqMat(t, k, pntcon, alpha, lagrange, E)
% [Deq veq] = BuildDeqMat(t, k, pntcon) or
% [Deq veq] = BuildDeqMat(t, k, pntcon, alpha, lagrange, E)
%
% Initialize shape equality constraint matrix "Deq" and rhs "veq" such that
% the spline point-wise constraints is formulated as
%      Deq*alpha = veq, (alpha beeing the vector of coeffs)

% Note: Deq is sparse, veq is full
% shape will contains modified shape with initialization
%
% [Deq veq Xeq] = BuildDeqMat(...) return a scaling diagonal matrix Xeq
% such that
% Xeq*Deq is normalized. The shape-constraints can be equivalently handled
%        (Xeq*Deq) * alpha = Xeq*veq
% for better numerical stability
%
% If alpha, lagrange, and E inputs are provided, the knot-derivative
% is returned and its adjoint
%   [Dt Dtt] = BuildDeqMat(t, k, pntcon, alpha, lagrange, E)
%

global BSFK_DISPLAY

% current number of basis splines
n = size(t,1) - k;

% number of shape constraints
l = length(pntcon);
P = [pntcon.p];
if any(P>=k)
    ibad = find(P>=k,1,'first');
    error(['BSFK: point-wise derivative order p(=%d) ' ...
           'must be smaller than k(=%d)'], P(ibad), k);
end

%% we will build the rhs of the equality constraints

if nargout>=2 && nargin<=3 % rhs is requested
    % just concatenate all the values together
    out2 = cat(1,pntcon.v);
end

%% Build the linear equality matrix D
% The result in this step does depend on the knot locations

% number of constraints
m = sum([pntcon.nc]);
    
% Building the equality matrix D
if nargin<=3
    
    % number of non-zeros elements of D
    % for each structure, there is m rows
    % each row has k non zeros elements (regardless p)
    nzD = m*k;
 
    % Loop on all constraint derivative
    % building the indexes and value of the linear constraint matrix
    ilast = 0;
    row = zeros(nzD,1,'double');
    col = zeros(nzD,1,'double');
    val = zeros(nzD,1,'double'); % sparse does not support single 
    xn = zeros(m,1,'double');
    rowstart = 0;
    singtol = eps(1)^2;
    for i=1:l
        s = pntcon(i);
        % Dp is sparse of the size mi x n
        Dp = Bspline(s.x, t, k, [], s.p);
        % it should find (n-p)*(p+1) elements
        [r c dp] = find(Dp);
        
        nz = length(dp); % == (n-p)*(p+1)
        mi = size(Dp,1); % == s.nc
        
        % add
        r = rowstart+r;
        idx = ilast+(1:nz);
        row(idx) = r;
        col(idx) = c;
        val(idx) = dp;
        dn = normest(Dp);
        if BSFK_DISPLAY && dn<singtol
            warning('BSFK:SingularCon', ...
                    'shape-constrained matrix close to singular');
        end
        xn(rowstart+(1:mi)) = 1/dn;
        % update the counters
        rowstart = rowstart+mi;
        ilast = ilast + nz;
    end % for-loop for D
    
    if any(isinf(xn))
        error('BSFK: singular shape-constrained matrix');
    end
    
   % Build the matrix for lower constraints
    if ilast<nzD
        idx = (1:ilast);
        Deq = sparse(row(idx),col(idx),val(idx),m,n);
    else % no need for subindexing
        Deq = sparse(row,col,val,m,n);
    end
   
    % Left diagonal scaling matrix (m x m)
    if nargout>=3
        Xeq = spdiags(xn, 0, m, m);
    end
else % Build the derivative of D with respect to knots
     % Build the knot-derivative of the matrix
    
    % Loop on all constraint derivative
    % building the indexes and value of the linear constraint matrix
    ne = size(E,2);
    rowstart = 0;
    DL = zeros([n ne], 'double');
    DR = zeros([m ne], 'double');
    
    lagrange = reshape(lagrange, [], 1);
    for i=1:l
        s = pntcon(i);
        p = s.p;
        mi = s.nc;
        idx = rowstart+(1:mi);
        u = lagrange(idx);
        % DDp is full of the size (n-p) x size(E,2)
        c = BsplineKnotDeriv(s.x, t, k, alpha, p, E, u);
        DR(idx,:) = c{1};
        % accumulate the gradient of the constraints
        DL = DL + c{2};
        rowstart = rowstart+mi;
    end % for-loop for D
    
    % m x size(E,2)
    % return the knot-derivative of the matrix as the first output
    Deq = DR;
    
    % n x size(E,2)
    % return the dual knot-derivative (Dtt) of the matrix as the second
    % output
    out2 = reshape(DL, [n ne]);
    
end % if

end % BuildDeqMat

%%
function [Deq out2 Xeq] = BuildPeriodicMat(t, k, alpha, lagrange, E)
% [Deq veq Xeq] = BuildPeriodicMat(t, k)
%
% Initialize shape equality constraint matrix "Deq" and rhs "veq" such that
% the spline periodic constraints is formulated as
%      Deq*alpha = veq, (alpha beeing the vector of coeffs)

% Note: Deq is sparse, veq is full
% shape will contains modified shape with initialization
%
% [Deq veq Xeq] = BuildDeqMat(...) return a scaling diagonal matrix Xeq
% such that
% Xeq*Deq is normalized. The shape-constraints can be equivalently handled
%        (Xeq*Deq) * alpha = Xeq*veq
% for better numerical stability
%
% If alpha, lagrange, and E inputs are provided, the knot-derivative
% is returned and its adjoint
%   [Dt Dtt] = BuildPeriodicMat(t, k, alpha, lagrange, E)
%

global BSFK_DISPLAY

% current number of basis splines
n = size(t,1) - k;

% number of constraints
l = k-1; % 0th, 1th..., (k-2)th derivatives
P = 0:k-2;

%% we will build the rhs of the equality constraints

if nargout>=2 && nargin<=2 % rhs is requested
    % derivative jumps are zeros
    out2 = zeros(l,1,'double');
end

%% Build the linear equality matrix D
% The result in this step does depend on the knot locations

% number of constraints
m = l;
    
% Building the equality matrix D
if nargin<=2
    
    % number of non-zeros elements of D
    % there is m rows
    % each row has respectively 2 x (1, 2, ..., m) non-zero elements
    nzD = m*(m+1);
 
    % Loop on all constraint derivative
    % building the indexes and value of the linear constraint matrix
    ilast = 0;
    row = zeros(nzD,1,'double');
    col = zeros(nzD,1,'double');
    val = zeros(nzD,1,'double'); % sparse does not support single 
    xn = zeros(m,1,'double');
    rowstart = 0;
    singtol = eps(1)^2;
    for p=P

        % Dp is sparse of the size 1 x n
        Dp = + Bspline(t(1), t, k, [], p, false) ...
             - Bspline(t(end), t, k, [], p, true);
        
        % it should find (n-p)*(p+1) elements
        [r c dp] = find(Dp);

        nz = length(dp); % == (n-p)*(p+1)
        mi = size(Dp,1); % == 1
        
        % add
        r = rowstart+r;
        idx = ilast+(1:nz);
        row(idx) = r;
        col(idx) = c;
        val(idx) = dp;
        
        dn = normest(Dp);
        if BSFK_DISPLAY && dn<singtol
            warning('BSFK:SingularCon', ...
                    'periodic-constrained matrix close to singular');
        end
        xn(rowstart+(1:mi)) = 1/dn;
        % update the counters
        rowstart = rowstart+mi;
        ilast = ilast + nz;
    end % for-loop for D
    
    if any(isinf(xn))
        error('BSFK: singular periodic-constrained matrix');
    end
    
   % Build the matrix for lower constraints
    if ilast<nzD
        idx = (1:ilast);
        Deq = sparse(row(idx),col(idx),val(idx),m,n);
    else % no need for subindexing
        Deq = sparse(row,col,val,m,n);
    end
   
    % Left diagonal scaling matrix (m x m)
    if nargout>=3
        Xeq = spdiags(xn, 0, m, m);
    end
else % Build the derivative of D with respect to knots
    % Build the knot-derivative of the matrix
    
    % Loop on all constraint derivative
    % building the indexes and value of the linear constraint matrix
    ne = size(E,2);
    rowstart = 0;
    DL = zeros([n ne], 'double');
    DR = zeros([m ne], 'double');
    
    lagrange = reshape(lagrange, [], 1);
    for p=P
        mi = 1;
        idx = rowstart+(1:mi);
        u = lagrange(idx);
        % DDp is full of the size (n-p) x size(E,2)
        % Dp is sparse of the size 1 x n
        cleft = BsplineKnotDeriv(t(1), t, k, alpha, p, E, u, false);
        cright = BsplineKnotDeriv(t(end), t, k, alpha, p, E, u, true);
        % cleft - cright
        c = cellfun(@minus,cleft,cright,'UniformOutput',false);
        DR(idx,:) = c{1};
        % accumulate the gradient of the constraints
        DL = DL + c{2};
        rowstart = rowstart+mi;
    end % for-loop for D
    
    % m x size(E,2)
    % return the knot-derivative of the matrix as the first output
    Deq = DR;
    
    % n x size(E,2)
    % return the dual knot-derivative (Dtt) of the matrix as the second
    % output
    out2 = reshape(DL, [n ne]);
    
end % if

end % BuildPeriodicMat

%%
function [C h] = BuildCMat(p, t)
% [C h] = BuildCMat(p, t)
% Build the constraint matrix C*t(p) >= h
% so that the knots do not bump into each other, this is also related
% to the zero-symmetrical box-constraint on Jupp's transformation:
%   -b <= kappa_i := log[ (t_i+1 - t_i) / (t_i-t-i-1) ] <= b

% Boor & Rice constant, must be small than 1
% if minratio 
minratio = 1/16; % 0.0625

% q will be unique and sorted
q = p;
q = union(p-1,q);
q = union(p+1,q);

loc = ismembc2(p, q);
l = length(p);
nz = 8*l;
col = zeros(1,nz,'double');
row = zeros(1,nz,'double');
val = zeros(1,nz,'double'); % sparse does not support single 
for i=1:l
    c = loc(i);
    % eqt 2.18 (a)
    i1 = (i-1)*8+(1:4);
    col(i1) = [c-1 c c-1 c+1];
    row(i1) = (1+(i-1)*2) + zeros(1,4);
    val(i1) = [-1 1 minratio -minratio];
    % eqt 2.18 (b)
    i2 = (i-1)*8+(5:8);
    col(i2) = [c+1 c c-1 c+1];
    row(i2) = (2+(i-1)*2) + zeros(1,4);
    val(i2) = [1 -1 minratio -minratio];
end

% knots constraints are: S*t(q) >= 0
S = sparse(row,col,val);

% Column of p
J = false(1,size(S,2));
J(loc) = true;

% knots constraints are; C*t(p) >= h
C = S(:,J);
tcp = -t(q(~J)); % minus fixed-knots
h = S(:,~J)*tcp(:);

end % BuildCMat

%%
function smoothing = UpdateConstraints(smoothing, t, ...
                                       shape, pntcon, periodic)
% smoothing = UpdateConstraints(smoothing, t, shape, pntcon, periodic)
% Update the constraint matrices and rhs when knots positions change
% Do not use after knot removal

    if nargin<3
        shape = smoothing.shape;
        pntcon = smoothing.pntcon;
        periodic = smoothing.periodic;
    end
    
    smoothing.t = t;
    
    knotidx = smoothing.knotidx;
    k = smoothing.k;
    
    [D LU X] = BuildDineqMat(t, knotidx, k, shape);
    
    % Put the shape-constrained information together smoothing
    smoothing.shape = shape;
    smoothing.D = D;
    smoothing.LU = LU;
    smoothing.X = X;
    
    [Deq veq Xeq] = BuildDeqMat(t, k, pntcon);
     
    if periodic
        [Deq_periodic veq_periodic Xeq_periodic] = BuildPeriodicMat(t, k);
        % Concatenate the point-wise constraints and periodic contraints
        % together
        Deq = cat(1, Deq, Deq_periodic);
        veq = cat(1, veq, veq_periodic);
        Xeq = blkdiag(Xeq, Xeq_periodic);
    end
    
    % Put the point-wise constrained information together the rest
    smoothing.pntcon = pntcon;
    smoothing.Deq = Deq;
    smoothing.veq = veq;
    smoothing.Xeq = Xeq;
   
end % UpdateConstraints

%%
function smoothing = InitPenalization(y, t, k, d, lambda, p, regmethod, ...
                                      knotidx, shape, pntcon, periodic)
% InitPenalization(y, t, k, d, lambda, p, regmethod, knotidx, ...
%                  shape, pntcon, periodic)
% Initialize smoothing and shape-constraints.
% To be called when change in knots topology occurs

% Initialize smoothing structure
smoothing = InitSmoothing(y, t, k, d, lambda, p, regmethod, knotidx, ...
                          periodic);
smoothing = UpdateConstraints(smoothing, t, shape, pntcon, periodic);
    
end % InitPenalization

%%
function [r J yfit OK] = BuildJacobian(x, t, k, GPflag, smoothing, ...
                                       shape, pntcon, periodic, dt)
% [r J] = BuildJacobian(x, t, k, GPflag, smoothing, shape, pntcon, ...
%                       periodic)
% 
% Build the Jacobian matrix (with respect to free knot) matrix for
% Gauss-Newton iteration
% r: is the fit residual
% J is the jacobian of r with respect to free-knots
% NOTES:
% [r g] = BuildJacobian(..., dt) compute the directional derivative
%   along the direction dt, i.e., g = Jacobian*dt
% [r J yfit] = BuildJacobian(...) return also fitting spline
% [r J yfit OK] = BuildJacobian(...), OK returned is false if Jacobian
% cannot be computes (colalision of overly constrained)

OK = false;

% Update the constraints matrices corresponding to new knot positions
smoothing = UpdateConstraints(smoothing, t, shape, pntcon, periodic);

if nargin<9 % direction of the derivative is not provided
    % Full Jacobian
    dt = smoothing.E;
end
                         
% Compute the basis and derivatives with respect to free-knot
% B (m x n): Model matrix
% dB (m x n x l): Knot-erivative of B
% Pseudo-inverse of B: Bi (n x m)
% alpha optimal coefficients (for the current set of knots) (n x 1)
% yfit: fitting data (m x 1)
% r: fit residual % (m x 1), r=P*ytarget with P = (speye(m)-B*Bi), a
%    projection matrix on orth span(B), note P is autoadjoint (symmetric)
[B Bi alpha yfit r dB] = ModelMat(x, t, k, smoothing, dt);

% coalsing knots!!
if isempty(B)
    [r J] = deal([]);
    return;
end

% m = number of data (including regularization)
% n = length(t)-k: number of basis function
% l = number of free knots
[m n l] = size(dB);

% 
% % Compute the derivative of the Moore-Penrose derivative
% % Golub & Oereyra differentiation formula
% Pm = B*Bi;
% Pn = Bi*B;
% Qm = eye(m)-Pm;
% Qn = eye(n)-Pn; % ~ zeros is spline matrix B is full rank
%                 % (i.e., Schoenberg-Whitney condition is satisfied)
% Bit = Bi.';
% dBt = dB.';
% dBi = -Bi*dB*Bi + Qn*(dBt*Bit)*Bi + (Bi*Bit)*dBt*Qm;

D = smoothing.D;
Deq = smoothing.Deq;
if isempty(D) && isempty(Deq)
    %% Unconstrained Jacobian with respect to free-knots
    % Kaufman
    a = reshape(alpha, [1 n]);
    dBa = bsxfun(@times, dB, a);
    dBa = reshape(sum(dBa,2), [m l]);
    % Kaufman matrix * y
    Ay = dBa - B*(Bi*dBa); % (m x l)
    
    if GPflag
        %% Transposed Kauffman
        rdB = bsxfun(@times, r, dB); % (m x n x l)
        rdB = reshape(sum(rdB,1), [n l]); % (n x l)
        Aty = (rdB.'*Bi).'; % (m x l)
        
        % Golub-Peyreyra
        J = -(Ay + Aty);
    else
        J = -Ay;
    end
else
    %% Shape-constrained Jacobian
    % Kaufman
    
    % Active set
    actset = Bi.actset;
    % Force the lagrange parameters to be zero on non active indices
    uineq = zeros(size(Bi.lagrange_ineq), class(Bi.lagrange_ineq));
    uineq(actset) = Bi.lagrange_ineq(actset);
    ueq = Bi.lagrange_eq;

    nact = size(actset,1)+size(Bi.lagrange_eq,1);
    if nact>=n
        J = nan(m,l);
        warning('BSFK:overconstrained', ...
                'BSFK: the spline seems to be overly constrained');
        return;
    end
    
    % Derivative of the active constraints wrt alpha 
    R = -[D(actset,:); ...
          Deq]; % nact x n
    
    % Derivative of the active constraints wrt knots 
    [Gamma_ineq Rtu_ineq] = BuildDineqMat(t, smoothing.knotidx, k, ...
                                          shape, alpha, uineq, dt);
    % number of pointwise constraints
    mpntcon = sum([pntcon.nc]);
    % The head of the dual parameter
    upntcon = ueq(1:mpntcon);
    % Derivative of the equality constraints wrt knots 
    [Gamma_eq Rtu_eq] = BuildDeqMat(t, k, pntcon, alpha, upntcon, dt);
    % Derivative of the periodic constraints wrt knots
    if periodic
        % The tail of the dual parameter
        iperiodic = mpntcon+1:size(ueq,1);
        uperiodic = ueq(iperiodic);
        [Gamma_periodic Rtu_periodic] = BuildPeriodicMat(t, k, alpha, ...
                                                   uperiodic, dt);
        % Concatenate the point-wise constraints and periodic contraints
        % together
        Gamma_eq = [Gamma_eq; ...
                    Gamma_periodic];
        Rtu_eq = Rtu_eq + Rtu_periodic;
    end
    Gamma = -[Gamma_ineq(actset,:); ...
              Gamma_eq]; % nact x l
    % Rtu: n x l
    Rtu = Rtu_ineq + Rtu_eq;
    
    % null space of R
    % - to try: scalled with X before for stability?
    N = qrnull(R); % n x (n-nact)
    % N = null(full(R)); % <- svd costly!
    
    BN = B*N; % m x (n-nact)
    BNi = pseudoinverse(BN); % (n-nact) x m
    
    if size(dt,2)>1
        a = reshape(-alpha, [1 n]);
        dBa = bsxfun(@times, dB, a); % (m x n x l)
        dBa = reshape(sum(dBa,2), [m l]); % (m x l)

        Ri = pseudoinverse(R); % (n x nact)
        dB2 = B*(Ri*Gamma); % (m x l)
        Psi = dBa + dB2; % (m x l)

        % Kaufman matrix
        Ay = Psi - BN*(BNi*Psi); % (m x l)

        if GPflag
            %% Transposed Kauffman
            rdB = bsxfun(@times, r, dB); % (m x n x l)
            rdB = reshape(sum(rdB,1), [n l]); % (n x l)
            K = rdB + Rtu; % (n x l)
            Phi = ((K.' * N) * BNi).'; % m x l

            % Golub-Peyreyra
            J = Ay - Phi;
        else
            J = Ay;
        end
    else % single directional derivative (typically when called from simul)
        dBa = dB*(-alpha); % (m x 1)

        Ri = pseudoinverse(R); % (n x nact)
        dB2 = B*(Ri*Gamma); % (m x 1)
        Psi = dBa + dB2; % (m x 1)

        % Kaufman matrix
        Ay = Psi - BN*(BNi*Psi); % (m x 1)
        
        if GPflag
            % Transposed Kauffman
            Kt = (r.'*dB) + Rtu.'; % (1 x n)
            Phi = ((Kt * N) * BNi).'; % (m x 1)

            % Golub-Peyreyra
            J = Ay - Phi;
        else
            J = Ay;
        end
    end
    
end % if on unconstrained and shape-constrained

% Signal Jacobian sucessully built
OK = true;

end % BuildJacobian


%%
function [s g] = GaussNewtonStep(t, p, r, J, C, h, qpengine)
%  [s g] = GaussNewtonStep(t, p, r, J, C, h)
% Solve the reduced problem by quadratic programming
% s: is increment of t (decent direction)
% g = J.'*r

global BSFK_DISPLAY
global BSFK_QPENGINE

% By default, use the global QP engine set by the main function
if nargin<7 || isempty(qpengine)
    qpengine = BSFK_QPENGINE;
end

% Hessian and residual gradient vector for QP
H = J.'*J;
% Symmetrize (eigs 'sa' need true symmetric matrix)
H = 0.5*(H.'+H); 
g = J'*r;

% Detect ill-conditioning
if size(H,1)>1
    % There could be a problem with EIGS as reported here
    % http://www.mathworks.fr/support/solutions/en/data/1-1HS9R3/index.html?product=ML&solution=1-1HS9R3
    try
        % For some reason EIGS throw time to time an error
        % "Error with ARPACK routine dneupd:
        % dnaupd did not find any eigenvalues to sufficient accuracy"
        eigsopt = struct('disp',0,'tol',1e-3);
        smallestev =  eigs(H,1,'sa',eigsopt);
    catch %#ok
        err = lasterror; %#ok
        % Check we bump to this annoying error 
        if ~isempty(strfind(err.message, 'ARPACK'))
            smallestev = 0; % must go on
        else
            rethrow(err);
        end
    end
    
    largestev = normest(H);
    if smallestev <= eps(largestev)
        % Add a small l2 regularization (small diagonal matrix) to make
        % the matrix better conditioned
        deflation = 2*max(-smallestev,eps(largestev));
        l = size(H,1);
        % H = H + diag(deflation*ones(l,1));
        H(1:l+1:end) = H(1:l+1:end) + deflation;
    end
end

% free knots
tf = reshape(t(p),[],1);
hr = h - C*tf; % new constraint rhs
eq = false(size(C,1),1);

% Call the QP routine to solve constrained quadratic problems
% Minimize J(s) := 1/2 s'*H*s + g'*s
%          such that C*s >= hr
switch lower(qpengine)
    case {1, 'quadprog'}, % mlint does not like mixed type
        % Blind coded by the author, who doesn't own optimization toolbox
        qpoptions = optimset('Display','off',...
                             'LargeScale','off');
        s0 = zeros(size(g));
        Aeq = []; beq = [];
        lb = []; ub = [];
        [s trash ier] = quadprog(H, g, -C, -hr, ...
                                 Aeq, beq, lb, ub, s0, qpoptions); %#ok
        if ier<0 && BSFK_DISPLAY>=1
            warning('BSFK:QPNonConvergence', ...
                'GaussNewtonStep: not convergence encountered by MINQDEF');
        end
    case {2, 'qpas' 'qpc'}, % mlint does not like mixed type
        qpoptions = struct('engine','qpas',...
                           'prt',false);
        [s dual ier] = qpmin(g, H, C, hr, eq, qpoptions); %#ok
        if ier && BSFK_DISPLAY>=1
            warning('BSFK:QPNonConvergence', ...
                'GaussNewtonStep: not convergence encountered by QPMIN');
        end
    otherwise % MINQDEF engine, written in Matlab
        qpoptions = {};
        prt = false;
        [s dual ier]=minqdef(g, H, C, hr, eq, prt, qpoptions{:}); %#ok
%         if ier && BSFK_DISPLAY>=1
%             warning('BSFK:QPNonConvergence', ...
%                 'GaussNewtonStep: not convergence encountered by MINQDEF');
%         end
end % switch

end % GaussNewtonStep


%%
function [f g indic yfit] = simul(indic, gamma, x, t0, p, s, k, ...
                                  smoothing, shape, pntcon, periodic)
% [f g indic yfit]=simul(indic, gamma, x, t0, p, s, k, smoothing, ...
%                         shape, pntcon)
% Simuation gateway for the linesearch nlis0

global BSFK_DISPLAY

% direction where the derivative will be computed
dt = expandfreeknot(s, t0, p);
t = t0 + gamma*dt;

% Update the constraints matrices corresponding to new knot positions
smoothing = UpdateConstraints(smoothing, t, shape, pntcon, periodic);
simuGPflag = true; % exact derivative

[r J yfit OK] = BuildJacobian(x, t, k, simuGPflag, smoothing, ...
                              shape, pntcon, periodic, dt);
% coalition of knots
if ~OK
    % save('BSFK_debug.mat','x','y','-mat');
    if BSFK_DISPLAY>=2
        fprintf('BSFK: unexpected coalesing knots\n');
    end
    % signal to nlis0 something is wrong
    indic = -1;
    f = NaN;
    g = NaN;
else

    % square of L2 residual
    f = 0.5*(r.'*r);

    % Gradient of f (derivative with respect to gamma)
    g = J.' * r;
end

end % simul

%%
function ps = prosca(x1, x2, varargin)
% function ps = prosca(x1, x2, varargin)
% Generic L2 scalar product
ps = x1.'*x2;
end % prosca

%%
function t=dcube(t, f, fp, ta, fa, fpa, tlower, tupper)
%
% function t=dcube(t, f, fp, ta, fa, fpa, tlower, tupper)
%
% Using f and fp at t and ta, computes new t by cubic formula
% safeguarded inside [tlower,tupper].
%
%  Original: Bruno Luong - 28/03/05
%

z1=fp+fpa-3*(fa-f)/(ta-t);
b=z1+fp;

%
% first compute the discriminant (without overflow)
%

if (abs(z1) <= 1)
   discri=z1*z1-fp*fpa;
else
   discri=fp/z1;
   discri=discri*fpa;
   discri=z1-discri;
   if ~xor(z1 >= 0, discri >= 0)
      discri=z1*discri;
   else
      discri=-1;
   end
end

if (discri < 0)
    if (fp < 0)
        t=tupper;
    else
        t=tlower;
    end
    %if (fp < 0) t=tupper; end
    %if (fp >= 0) t=tlower; end
    %t=min(max(t,tlower),tupper);
    return
end

%
%  discriminant nonnegative, compute solution (without overflow)
%

discri=sqrt(discri);

if (t-ta < 0)
    discri=-discri;
end

sign=(t-ta)/abs(t-ta);
if (b*sign > 0)
  t=t+fp*(ta-t)/(b+discri);
else
  den=z1+b+fpa;
  anum=b-discri;
  if (abs((t-ta)*anum) < (tupper-tlower)*abs(den))
    t = t+anum*(ta-t)/den;
  else
    t = tupper;
  end
end

t=min(max(t,tlower),tupper);

end % dcube

%%
function [t xn fn g logic nap] = ...
         nlis0(n, simul, prosca, xn, fn, fpn, t, tmin, tmax, d, g, ...
               amd, amf, imp, io, nap, napmax, varargin) %#ok
%
% function [t, xn fn g logic nap] = ...
%         nlis0(n, simul, prosca, xn, fn, fpn, t, tmin, tmax, d, g, ...
%               amd, amf, imp, io, nap, napmax, varargin)
%
%  linear search for minimization of h(t):=f(x + t.d), t>=0;
%  cubic interpolation + Wolfe's tests
%
%  INPUT:
%     n: dimension of the full space
%     simul: simulator function handle
%     prosca: scaler product function handle
%     xn(n): current vector
%     fn: current cost function
%     fpn: derivative of h(t) at t=0
%     t: first guess of t
%     tmin, tmax: t-bracket interval
%     d(n): direction of linear search
%     g(n): current gradient vector (of f at xn)
%     amd, amf: constants for Wolfe's tests
%     imp, io: level and of file identifier for output printing
%     nap: current number of simulator calling
%     napmax: maximum number of allowed simulator calling
%
%  OUTPUT:
%     t: new t-value
%     xn(n): new current vector
%     fn: new function value
%     g(n): new gradient
%     logic =
%        0          line search OK
%        1          minimization stucks
%        4          nap > napmax
%        5          user stopping
%        6          function and gradient do not agree
%        < 0        contrainte implicite active
%
%     nap: updated number of simulator calling
%
%  Original: Bruno Luong - 28/03/05
%            29/10/05: output messages according to imp
%
     
if ~(n>0 && fpn<0 && t>0 && tmax>0 && amf>0 ...
     && amd>amf && amd<1)
  logic=6;
  return
end

tesf=amf*fpn;
tesd=amd*fpn;
barmin=0.01;
barmul=3;
barmax=0.3;
barr=barmin;
td=0;
tg=0;
fg=fn;
fpg=fpn;
ta=0;
fa=fn;
fpa=fpn;

ps=prosca(d,d,varargin{:});
d2=ps;

t=max(t,tmin);
if (t > tmax)
    tmin=tmax;
    if imp>=0
        fprintf('tmin forced to tmax\n');
    end
end

while (fn+t*fpn >= fn+0.9*t*fpn)
  t=2*t;
end

if imp>=4
  fprintf('  lsearch  fpn=%1.3e d2=%1.2e tmin=%1.2e tmax=%1.2e\n', ...
          fpn, d2, tmin, tmax);
end

%
% Label 30
%
indica=1;
logic=0;
if (t > tmax)
    t=tmax;
    logic=1;
end

x=xn+t*d;

% Main Loop, label 100
while 1
    
  nap=nap+1;
  if (nap > napmax)
    logic=4;
    fn=fg;
    xn=xn+tg*d;
    return
  end
  
  indic=4;
  %
  % Call the simulator
  %
  [f,g,indic]=simul(indic,x,varargin{:});
  if (indic == 0) % User stopping
    logic=5;
    fn=f;
    xn=x;
    return;
  elseif (indic < 0) % Simulator can't perform
    td=t;
    indicd=indic;
    logic=0;
    if imp>=4
      fprintf(['      lsearch: simulator can''t perform at '...
               't=%1.3e indic=%d\n'], t, indic);
    end
    t=tg+0.1*(td-tg);
    %go to 905
  elseif (indic > 0) % Simulator OK, let's go...
      
    ps=prosca(d,g,varargin{:});
    fp=ps;
  %
  % First Wolfe's test
  %
    interpflag=0; %#ok
    ffn=f-fn;
    if (ffn > t*tesf)
      td=t;
      fd=f;
      fpd=fp;
      indicd=indic;
      logic=0;
      interpflag=1;
      if imp>=4
        fprintf('    lsearch: %1.3e %1.3e %1.3e\n', t, ffn, fp);
      end
      % go to 500
    else
  %
  % First test OK, perform second Wolfe's test
  %
      if imp>=4
        fprintf('    lsearch: %1.3e %1.3e %1.3e\n', t, ffn, fp);
      end
      if (fp > tesd) || (logic ~= 0)
        if (fp > tesd)
            logic=0;
        end
        fn=f;
        xn=x;
        return;
      end
  
  %
  % Label 350
  %
      tg=t;
      fg=f;
      fpg=fp;
      interpflag=(td ~= 0);
    end % if (ffn > t*tesf)
  
%
% extrapolation
%
    if ~interpflag
      taa=t;
      gauche=(1+barmin)*t;
      droite=10*t;
      t=dcube(t, f, fp, ta, fa, fpa, gauche, droite);
      ta=taa;
      if (t >= tmax)
        logic=1;
        t=tmax;
      end
    else % interpflag==1
%
% Label 500, interpolation
%
      if (indica <= 0)
        ta=t;
        t=0.9*tg+0.1*td;
      else
        test=barr*(td-tg);
        gauche=tg+test;
        droite=td-test;
        taa=t;
        t=dcube(t, f, fp, ta, fa, fpa, gauche, droite);
        ta=taa;      
        if (t > gauche && t < droite)
          barr=barmin;
        else
          barr=min(barmul*barr,barmax);
        end
      end
    end % if ~interpflag

%
% --- fin de boucle
%     - t peut etre bloque sur tmax
%       (venant de l'extrapolation avec logic=1)
%
% label 900
%
    fa=f;
    fpa=fp;
    
  end % if (indic == 0 ...)
%
% label 905  
%
  indica=indic;
%
% do we have to pursue ?
%

  if (td == 0)
      exitloop=0;
  elseif (td-tg < tmin)
      exitloop=1;
  else
%
%     --- limite de precision machine (arret de secours) ?
%
      z=xn+t*d;
      exitloop=all(z==xn | z==x);
  end
  
  if exitloop
    logic=6;
%
%     si indicd<0, derniers calculs non faits par simul
%
    if (indicd < 0)
      logic=indicd;
    end
%
%     si tg=0, xn = xn_depart,
%     sinon on prend xn=x_gauche qui fait decroitre f
%
    if (tg ~= 0)
      fn=fg;
      xn=xn+tg*d;
    end
% label  940
    if imp>0
      fprintf('    lsearch: %1.3e %1.3e %1.3e\n', tg, fg, fpg);
      if logic==6
        fprintf('    lsearch: %1.3e %1.3e %1.3e\n', td, fd, fpd);
      elseif logic==7
        fprintf('    lsearch: %1.3e indic=%d\n', td, indicd);  
      end
    end
    return
  end
%
% recopiage de x et boucle
%
% label 950
%
  x=xn+t*d;
  
end

end % nlis0

%%
function [gammaopt g0] = linesearch(x, t0, p, s, k, smoothing, ...
                                    shape, pntcon, periodic)
% [gammaopt g0] = linesearch(x, t0, p, s, k, smoothing, shape, pntcon, periodic)
% gammaopt is the "optimal" step such that t0 + gammaopt*s decreases
% in a "good" way

% Adjustable constant, we recommend 1/8 - larger than value used in
% n1qn3
delta = 1/8;
%
% Constants used in two Wolfe's tests
%
rm1=delta;
rm2=0.9;

l = length(p);

imp = 0;
io = 1;
nap = 1;
napmax = 10*l;

gmin = 0;
gmax = 1;
gopt = 15/16; % 0.9375

simparams = {x, t0, p, s, k, smoothing shape pntcon periodic};
gamma = 0;
indic = 1;
[f0 g0] = simul(indic, gamma, simparams{:});

n = 1;
d = 1;
fpn=prosca(d,g0,simparams{:});
[gammaopt trash fopt gopt logic] = ...
    nlis0(n, @simul, @prosca, gamma, f0, fpn, gopt, gmin, gmax, ...
               d, g0, rm2, rm1, imp, io, nap, napmax, ...
               simparams{:}); %#ok
           
 if ~(fopt<f0) %|| logic>0
     gammaopt = 0;
 end

end % linesearch
