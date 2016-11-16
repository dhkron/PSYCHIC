function varargout = Bspline(x, knots, k, alpha, rtab, varargin)
% [s1 s2 ...] = Bspline(x, knots, k, alpha, rtab)
%
% Compute the Bspline function with coefficient alpha
%
% INPUTS:
%   x: vectors/array, point abscissa coordinates at which the function is
%      to be evaluated
%   knot: vector, knots points, must be ascending sorted
%   k-1: is the order of the spline
%   alpha: coefficient of the B-splines
% OUTPUT:
%   s1, s2... : spline functions values/derivative computed at x
%       each si to rtab(i)
%
% when rtab = [r1, r2 ...] is provided, the rth derivatives (with respect
% to abscissa) are returned
%   [sr1 sr2...] = Bspline(x, knots, k, alpha, rtab)
% sr1(x) is d^(r1)s/dx^(r1) (x), etc...
%
%  Note: si has support on [knots(j),knots(j+k-r)[
%           r := rtab(i)
%
%  Call with TRUE for 6th argument: Bspline(..., TRUE) to compute
%  the "Left B-spline", which has support on ]knots(j),knots(j+k-r)].
%  The left B-spline is needed e.g. to evaluate recursion integral between
%  two B-splines.
% 

if nargin<4
    alpha = [];
end

if nargin<5 || isempty(rtab)
    rtab = 0;
end

x = x(:);
rtab = reshape(rtab, 1, []);

out = cell(1,size(rtab,2));
for i=1:length(rtab)
    r = rtab(i);
    if r>0
        [ap td kd] = DerivB(knots, k, r, alpha);
        s = Bernstein(x, td, [], kd, ap, varargin{:});
    else
        s = Bernstein(x, knots, [], k, alpha, varargin{:});
    end
    out{i} = s;
end

varargout = out;