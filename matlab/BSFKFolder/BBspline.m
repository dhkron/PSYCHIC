function varargout = BBspline(varargin)
% [s1 s2 ...] = BBspline(x, knots, k, alpha, rtab)
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
% Author: Bruno Luong <brunoluong@yahoo.com>
% History: 18-Nov-2009: Original
%

out = cell(1,nargout);
[out{:}] = Bspline(varargin{:}); % private method, because having the same
                                 % name with Matlab function
varargout = out;