function varargout = BsplineKnotDeriv(x, knots, k, alpha, rtab, dknots, ...
                                      lambda, varargin)
% c = BsplineKnotDeriv(x, knots, k, alpha)
%
% c = BsplineKnotDeriv(x, knots, k, alpha, rtab, dknots)
%
% c = BsplineKnotDeriv(x, knots, k, alpha, rtab, dknots, lambda)
%
% Compute the B-spline function having coefficient alpha and its derivative
% with respect to knots
%
% INPUTS:
%   x: vectors/array, point abscissa coordinates at which the function is
%      to be evaluated
%   knots: vector, knots points, must be ascending sorted
%   k-1: is the order of the spline
%   alpha: coefficient of the B-splines
%   dknots: directions of the knot derivative, ranged in p columns
%   lambda: dual variable
% OUTPUT:
%   c: is a 1x2 cell, where
%      c{1} spline function value computed at x
%      c{2} knot-derivative value computed at x in the direction dknots
%   if the dual variable lambda is provided
%      c{1} knot-derivative value computed at x in the direction dknots
%           right-product with alpha
%      c{2} knot-derivative value computed at x in the direction dknots
%           left-product with lambda
%
% when rtab = [r1, r2 ...] is provided, the rth derivatives (with respect
% to abscissa) are returned
%   [c1 c2...] = BsplineKnotDeriv(x, knots, k, alpha, rtab)
%     c1{1} is the spatial derivative d^(r1)s/dx^(r1) (x),
%     c1{2} is the knot-derivative of c1{1}
%   Similarly c2, c3, ... correspond respectively to the order r2, r3, ...
%
%  Call with TRUE for 8th argument: B = BsplineKnotDeriv(..., TRUE) to
%  compute the "Left B-spline", which has support on ]t(j),t(j+k)].
%
%  Last update: 31-May-2010: Add leftflag
% 

if nargin<5 || isempty(rtab)
    rtab = 0;
end

if isequal(size(dknots),[1 length(knots)])
    dknots = dknots.';
end

% reshape in column
x = x(:);

rtab = reshape(rtab, 1, []);
out = cell(1,numel(rtab));
if nargin<7
    for i=1:length(rtab)
        r = rtab(i);
        if r>0
            % ap: coefficients of the rth derivative
            % dap: knot-derivative of those coefficients
            [ap td kd dap subs] = DerivBKnotDeriv(knots, k, r, dknots, alpha);
            
            % knot-derivative of the basis
            dkd = dknots(subs,:);
            [sr dsr] = BernKnotDeriv(x, td, [], kd, dkd, ap);
            
            % Add the contribution of the spline built from the
            % knot-derivative of the coefficients
            dsr = dsr + Bernstein(x, td, [], kd, dap);
            
        else
            [sr dsr] = BernKnotDeriv(x, knots, [], k, dknots, alpha);
        end
        out{i} = {sr dsr};
    end
else % dual variable lambda is provided
    for i=1:length(rtab)
        r = rtab(i);
        if r>0
           
            subs = r+1:length(knots)-r;
            td = knots(subs);
            kd = k-r;
            
            % knot-derivative of the basis
            dkd = dknots(subs,:);
            [B dB] = BernKnotDeriv(x, td, [], kd, dkd, ...
                                   [], [], varargin{:});   
            
            % ap: coefficients of the rth derivative
            % dap: knot-derivative of those coefficients
            [Dr td kd dDr] = DerivBKnotDeriv(knots, k, r, dknots); %#ok
            
            % The derivarive of B*Dr            
            [m n p] = size(dB);
            n = n+r;
            
            % Compute left and right products: 
            %   ds = (B*dDr + dB*Dr)*alpha
            %   dalpha = lambda.'*(B*dDr + dB*Dr)
            lambda = lambda(:).'; % row-vector
            % These are simple matrix x vector
            lbdB = lambda*B; % 1 x (n-r)
            Dralpha = Dr*alpha; % (n-r) x 1
            
            % This can be carried out outside the for-lpp
            lbddB = lambda*reshape(dB,m,[]); 
            lbddB = reshape(lbddB, n-r, p); % (n-r) x p
            % Preallocate
             cls = class(B);
            % Preallocate
            dDralpha = zeros([n-r p],cls); % (n-r) x p
            lbdBdDr = zeros([n p],cls); % (n x p)
            dBDralpha = zeros([m p],cls); % (m x p)
            if ~iscell(dDr)
                dDr = {dDr};
            end
            for j=1:p
                dDralpha(:,j) = dDr{j}*alpha; % (n-r) x 1
                lbdBdDr(:,j) = lbdB*dDr{j}; % (1 x n)
                dBDralpha(:,j) = dB(:,:,j)*Dralpha; % (m x 1) 
            end
            % left product with the dual
            ds = B*dDralpha + dBDralpha; % (m x p)
            % right product with the coefficient
            dalpha = (Dr.')*lbddB + lbdBdDr; % (n x p)
            
%             dBD = zeros([m n p], cls);
%             for j=1:p
%                 dBD(:,:,j) = dB(:,:,j)*Dr + B*dDr{j};
%             end    
%             
%             % left product with the dual
%             lambda = lambda(:).'; % row-vector
%             dalpha = lambda*reshape(dBD, [m n*p]);
%             dalpha = reshape(dalpha, [n p]);
%             
%             % right product with the coefficient
%             alpha = alpha(:); % column-vector
%             % permute the coefficient dimension to the end
%             dBD = reshape(permute(dBD, [1 3 2]), [m*p n]);
%             ds = dBD*alpha;
%             ds = reshape(ds, [m p]);
            
        else
            [junk ds dalpha] = BernKnotDeriv(x, knots, [], k, dknots, ...
                                             alpha, lambda, varargin{:}); %#ok
        end
        out{i} = {ds dalpha};
    end
end

varargout = out;