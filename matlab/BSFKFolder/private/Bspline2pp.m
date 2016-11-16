function pp = Bspline2pp(sresult)
% function pp = Bspline2pp(sresult)
% Convert the knot-vector of B-spline to pp form
% History:
%   07-Jun-2010, scale the abscissa to make pp determination
%   more robust

% if nargin<2
%     extrapflag = false;
% end

k = sresult.k;
f = sresult.unnormfun;
t = reshape(sresult.t, 1, []);

knots = unique(t);
left = knots(1:end-1);
right = knots(2:end);
dt = right-left;

% Compute k values of B-spline on each sub-interval
a = linspace(0,pi,k+2).';
a = 0.5*(1-cos(a(2:end-1)));
xi = bsxfun(@times, a, dt);
x = bsxfun(@plus, left, xi);
y = Bernstein(x, t, [], k, sresult.alpha);
if k==1
    y = y.';
end

% Regression to find the polynomial coefficients
p = k-1:-1:0;
xi = repmat(a, [1 1 size(x,2)]);
xik = bsxfun(@power, xi, p);
coefs = MultiSolver(xik,y);
dt = f(right)-f(left);
s = bsxfun(@power, 1./dt, p(:));
coefs = s.*coefs;

% Dead code
% xi = bsxfun(@minus, f(x), f(left));
% xi = reshape(xi, [k 1 size(xi,2)]);
% xik = bsxfun(@power, xi, k-1:-1:0);
% coefs = MultiSolver(xik,y);

% Make dure breaks cover exactly the data
if knots(end-1) < 1
    knots(end)=1;
end
breaks = f(knots);

% Create Matlab pp-form
pp = struct('form', 'pp', ...
            'breaks', breaks, ...
            'coefs', coefs.', ...
            'pieces', length(left),...
            'order',k,...
            'dim',1);
            
end % Bspline2pp