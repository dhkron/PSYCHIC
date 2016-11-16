function varargout = testBSFK(caseid)
% function result = testBSFK(caseid)
%
% test BSFK package
% caseid: 1 -> function with three slopes, one shaped constrained
%              point-wise constraints on derivative and spline function
% caseid: 2 -> gaussian, two shaped constrained
%              one point-wise constraint on spline function
% caseid: 3 -> random, periodic
% otherwise -> random, unconstrained
%
% see also: BSFK, ppval
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% History: 17-Nov-2009: Original
%          19-Nov-2009: return result
%          04-Dec-2009: constrained test
%          08-Dec-2009: point-wise test
%          31-May-2010: Periodic test

if nargin<1
    caseid=0;
end

if isstruct(caseid)
    inputstruct = caseid;
    caseid = 100;
end

%%
switch caseid
    case 1,
        
        stdnoise = 0.02;
        
        x = linspace(0,3);
        y = x-1;
        y(x<=2) = 1;
        y(x<=1) = x(x<=1);
        
        y = y+stdnoise*randn(size(y));
        
        nknots = 20;
        fixknots = [];
        k = 4;
        clear shape

        % function increasing
        lo = 0;
        up = +inf;
        shape(1) = struct('p', 1, 'lo', lo, 'up', up);
        
        % Function forced to have follownh values
        %    s(0)=1, s(3/2)=1, s(2)=2
        %    s'(0)=1 s'(3)=1
        pntcon(1) = struct('p', 0, 'x', [0 1.5 3], 'v', [0 1 2]);
        pntcon(2) = struct('p', 1, 'x', [0 3], 'v', 1);
        
        options = struct('animation', 1, ...
            'figure', 1, ...
            'waitbar', 1, ...
            'display', 1, ...
            'd', 1, 'lambda', 0e-3, 'regmethod', 'c', ...
            'qpengine', '', ...
            'sigma', [], ...
            'shape', shape, ...
            'pntcon', pntcon);
        
    case 2,
        
        stdnoise = 0.02;
        
        x = linspace(-1,1);
        y = exp(-6*(x.^2));
        
        y = y+stdnoise*randn(size(y));

        nknots = [10 4 10];
        fixknots = [-0.18 0.18];
        k = 4;
        % function increasing
        lo = 0;
        up = +inf;
        shape(1) = struct('p', 0, 'lo', lo, 'up', up);
        % constraints on second derivative -> convex for |x|>=2
        lo = -inf(1,sum(nknots)); lo(1:10)=0; lo(end+(-9:0))=0;
        up = +inf;
        shape(2) = struct('p', 2, 'lo', lo, 'up', up);
        
        % zero-derivative of the central point x=0
        pntcon(1) = struct('p', 1, 'x', 0, 'v', 0);
        options = struct('animation', 1, ...
            'figure', 1, ...
            'waitbar', 1, ...
            'display', 1, ...
            'd', 2, 'lambda', 0e-6, 'regmethod', 'c', ...
            'qpengine', '', ...
            'sigma', [], ...
            'shape', shape, ...
            'pntcon', pntcon);
        
    case 3, %% Periodic
        
        %%
        stdnoise = 0.05;
        nknots=10;
        t=cumsum(rand(1,nknots+1));
        
        x=linspace(min(t),max(t)-eps(max(t)),2^9+1);
        
        %%
        k = 4; % C2
        lt = size(t,2);
        idx = [1+zeros(1,k) 2:lt-1 lt+zeros(1,k)];
        t = t(idx);
        % Number of basis
        n = length(t)-k;
        B = BBspline(x(:), t, k);
        
        alphatrue=randn(n,1);
        y = B*alphatrue;
        y = 1e3*(y+stdnoise*randn(size(y)));
        
        nknots = 20;
        fixknots = [];
        k = 4;
        options = struct('animation', 1, ...
            'figure', 1, ...
            'waitbar', 1, ...
            'display', 1, ...
            'd', 2, 'lambda', 0e-7, 'regmethod', 'c', ...
            'qpengine', '', ...
            'sigma', [], ...
            'periodic', true); 
    case 100,
        
        [x, y, k, nknots, fixknots, options] = deal(inputstruct.x, ...
            inputstruct.y, inputstruct.k, inputstruct.nknots, ...
            inputstruct.fixknots, inputstruct.options);
        
    otherwise
        %%
        stdnoise = 0.05;
        nknots=10;
        t=cumsum(rand(1,nknots+1));
        
        x=linspace(min(t),max(t)-eps(max(t)),2^9+1);
        
        %%
        k = 4; % C2
        lt = size(t,2);
        idx = [1+zeros(1,k) 2:lt-1 lt+zeros(1,k)];
        t = t(idx);
        % Number of basis
        n = length(t)-k;
        B = BBspline(x(:), t, k);
        
        alphatrue=randn(n,1);
        y = B*alphatrue;
        y = 1e3*(y+stdnoise*randn(size(y)));
        
        nknots = 20;
        fixknots = [];
        k = 4;
        options = struct('animation', 1, ...
            'figure', 1, ...
            'waitbar', 1, ...
            'display', 1, ...
            'd', 2, 'lambda', 0e-7, 'regmethod', 'c', ...
            'qpengine', '', ...
            'sigma', []);
       
end

pp = BSFK(x, y, k, nknots, fixknots, options);

fig = figure(2);
clf(fig);
ax = axes('Parent',fig);
minx = min(x);
maxx = max(x);
dx = (maxx-minx);
periodic = isfield(options, 'periodic') && options.periodic;
if periodic
    margin = 0.25;
else
    margin = 0;
end
left = minx - margin*dx;
right = maxx + margin*dx;
xi = linspace(left,right,1025);
if periodic
    xwrapped = minx + mod(xi-minx, dx);
    yfit =  ppval(pp,xwrapped);
else
    yfit =  ppval(pp,xi);
end
plot(ax, x, y, 'g.', ...
     xi, yfit,'r');

result = struct( ...
    'x', x, ...
    'y', y, ...
    'k', k, ...
    'nknots', nknots, ...
    'fixknots', fixknots, ...
    'options', options, ...
    'pp', pp);


if nargout>=1
    varargout = {result};
end
