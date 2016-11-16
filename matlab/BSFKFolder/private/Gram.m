function [R dR dt] = Gram(t, k, dt, freeidx)
% [R dR dt] = Gram(t, k, freeidx) OR
% [R dR dt] = Gram(t, k, [], dt)
%
% Compute the L2 Gram's matrix for B-spline basis of knots t
%
% INPUTS:
%   t: array of knots, must be ascending sorted (not checked)
%      elements could be repeated
%   k-1: the order of the spline
%   freeidx: index of (free) knots where the derivative of the Gram matrix
%            will be computed
%   dt: matrix of derivative directions (freeidx argument will be ignored)
%       size(dt,1) = size(t,1) (number of knots)
%       size(dt,2) = size(dR,3)
% OUTPUTS:
%   R: the Gram's matrix of basis, i.e.,
%   R(i,j) = Integral_dx * Bernstein(x,t,i,k)*Bernstein(x,t,j,k)
%   dR: the matrix of derivative; arranged by columns of dt
%       size(dR,3) = size(dt,2)
%       size(dR(:,:,i)) = size(R) for all i
%
% Implementation of the integral between two spline by recusion 
% method (see Boor/Lyche/Schumaker - On calculating with B-spline II.
% - Integrationn Numerische Methoden der Approximationtheorie, Vol. 3
% ISNM 30; or
% - Lyche/Morken/Quak, Theory and algorithms for non-uniform spline
% wavelets (lmq.pdf).
%
% RESTRICTION: knots t must have duplicated points in both ends so that the
%              formula of derivative with respect to abscissa coordinates
%              used here can be correctlty applied. Otherwise the result
%              will be wrong!!!!
% Author: Bruno Luong
%   10-Jun-2010: Change in Bernstein.m to avoid NaN for end knots

t = reshape(t,[],1);
nknots = size(t,1);
cls = class(t);

if k<1 || k>=nknots
    error('GRAM: k(=%d) must be positive and smaller than %d', k, nknots)
end

if nargin >=4
    % Default direction from freeidx
    r = max(k-1,0);
    if isempty(freeidx)
        freeidx = 1+r:nknots-r;
    else
        if any(freeidx<1+r | freeidx>nknots-r)
            error('GRAM: freeidx must comprise in (%d,%d)', 1+r, nknots-r);
        end
    end
    % vector of direction for derivative
    dt = zeros(nknots, numel(freeidx), cls);
    dt(freeidx + (0:numel(freeidx)-1)*nknots) = 1;
elseif nargin<3 || isempty(dt)
    dt = zeros(nknots,0);
elseif size(dt,1)~=nknots
    error('GRAM: dt must have %d rows', nknots);
end

[R dR] = quadengine(t, k, dt);

i = 1:nknots-k;
n = size(i,2);
p = size(dt,2);

% c = 1 / sqrt(nchoosek(2*k-1,k)*k);
% dt = c*dt;
% t = c*t;

tdiff = t(i+k)-t(i); % n x 1
dtdiff = reshape(dt(i+k,:)-dt(i,:), [n 1 p]); % n x 1 x p
ddtdifft = reshape(dtdiff, [1 n p]); % 1 x n x p

dtij = bsxfun(@times, tdiff, tdiff.'); % n x n
ddtij = bsxfun(@times, tdiff, ddtdifft) +  ...
        bsxfun(@times, dtdiff, tdiff.'); % n x n x p

dR1 = bsxfun(@times, ddtij, R); % n x n x p
dR2 = bsxfun(@times, dtij, dR); % n x n x p
% dR2(dtij==0) = 0;
dR = dR1 + dR2; % n x n x p
R = dtij.*R;

c = nchoosek(2*k-1,k)*k;
R = R / c;
dR = dR / c;

end % Gram


%%
function [R dR] = quadengine(t, r, dt)
% function [R dR] = quadengine(t, r, dt)

cls = class(t);

% k, l in [1,2...r]
[L K] = ndgrid(0:r,0:r); % <- [K(:) L(:)] is sorted by rows

% Generate array of nchoosek(K+L-1,K)
C=zeros(r+1,r+1);
C(2:end,1) = 1;
for j=1:r
    C(:,j+1) = cumsum(C(:,j));
end

nknots = size(t,1);
nkl = numel(K);

% This array will be used to access to third index of T, which contain
% combination of both k and l
idxkl = accumarray([K(:) L(:)]+1, (1:nkl).', [r+1 r+1]);
% (i,j,kl)
T = zeros([nknots nknots nkl],class(t));

p = size(dt,2);
% dT(:,:,:,i) is the derivative of T with respect to the direction dt(:,i)
dT = zeros([nknots nknots nkl p],cls);

% we'll need this array several times, so create it once
allknotidx = 1:nknots;
            
% Let's go!!!
for n = 2:nkl % for n=1, T should be filled with 0, so we don't do
    k = K(n);
    l = L(n);
    
    %%
    if k==0 % && l>0
        j = 1:nknots-l;
        
        % Derivative with respect to the second t
        % B: nknots x (nknots-l)
        % dB: nknots x (nknots-l) x p
        [B dB] = Btt(t, allknotidx, j, l, dt, true);
%         if any(isnan(dB(:)))
%             keyboard
%         end
        dtj = t(j+l)-t(j); % (nknots-l) x 1
        ddtj = dt(j+l,:)-dt(j,:); % (nknots-l) x p
        c = 1./dtj; % (nknots-l) x 1
        dc = bsxfun(@times,-c./dtj,ddtj); % (nknots-l) x p
        B = B.'; % (nknots-l) x nknots
        T(j,:,n) = bsxfun(@times, c, B); % (nknots-l) x nknots
        dc = reshape(dc,[nknots-l 1 1 p]);
        dB = reshape(permute(dB,[2 1 3]),[nknots-l nknots 1 p]);
        dT(j,:,n,:) = bsxfun(@times, dc, B) + ...
                      bsxfun(@times, c, dB); % (nknots-l) x nknots x 1 x p
        continue
    end
    %%
    if l==0 % && k>0
        i = 1:nknots-k;
        % B: nknots x (nknots-k)
        % dB: nknots x (nknots-k) x p
        [B dB] = Btt(t, allknotidx, i, k, dt, false);
        dti = t(i+k)-t(i); % (nknots-k) x 1
        ddti = dt(i+k,:)-dt(i,:); % (nknots-k) x p
        c = 1./dti; % (nknots-k) x 1
        dc = bsxfun(@times,-c./dti,ddti); % (nknots-k) x p
        c = c.'; % 1 x (nknots-k)
        T(:,i,n) = bsxfun(@times, c, B);  % nknots x (nknots-k)
        dc = reshape(dc,[1 nknots-k 1 p]);
        dB = reshape(dB,[nknots nknots-k 1 p]);
        dT(:,i,n,:) = bsxfun(@times, dc, B) + ...
                      bsxfun(@times, c, dB); % nknots x (nknots-k) x 1 x p
        continue
    end

    ck = C(n); % nchoosek(k+l-1,k);
    cl = (ck*k)/l; % nchoosek(k+l-1,l);
    
    jvec = (1:nknots-l).'; % (nknots-l) x 1
    %%
    for i=1:nknots-k
        dti = t(i+k)-t(i);       
        if dti==0
            j = jvec;
            dtj = t(j+l)-t(j); % (nknots-l) x 1
            idx = dtj~=0;
            if any(idx) % m = sum(idx)
                jidx = j(idx);
                m = size(jidx,1);
                dtjidx = dtj(idx); % m x 1
                ddtjidx = dt(jidx+l,:)-dt(jidx,:); % m x p
                c = ck./dtjidx; % m x 1
                dc = bsxfun(@times,-c./dtjidx,ddtjidx); % m x p
                % B: 1 x m
                % dB: 1 x m x p
                [B dB] = Btt(t, i, j(idx), l, dt, true);
                B = B.'; % m x 1
                T(jidx,i,n) = c.*B; % m x 1
                dc = reshape(dc,[m 1 1 p]);
                dB = reshape(dB,[m 1 1 p]);
                dT(jidx,i,n,:) = bsxfun(@times, dc, B) + ...
                                 bsxfun(@times, c, dB); % m x 1 x 1 x p
            end
        else % dti~=0
            
            dtj = t(jvec+l)-t(jvec); % (nknots-l) x 1
            idx = dtj==0;
            if any(idx) % m = sum(idx)
                jidx = jvec(idx);
                m = size(jidx,1);
                ddti = dt(i+k,:)-dt(i,:); % 1 x p 
                c = cl/dti; % 1 x 1
                dc =  (-c./dti).*ddti; % 1 x p
                % B: m x 1
                % dB: m x 1 x p
                [B dB] = Btt(t, jidx, i, k, dt, false);
                T(jidx,i,n) = c*B; % m x 1
                dc = reshape(dc,[1 1 1 p]);
                dB = reshape(dB,[m 1 1 p]);
                dT(jidx,i,n,:) = bsxfun(@times, dc, B) + ...
                                 c*dB; % m x 1 x 1 x p
            end
            
            loopidx = ~(idx);
            % Loop (on j) only when the two supports intersect 
            loopidx(i+k:end) = false;
            loopidx(1:i-l) = false;
            
            if k<=l
                
                iklm1 = idxkl(k+1,l);
                ikm1l = idxkl(k,l+1);
                Tvec1 = T(jvec,i,iklm1);
                Tvec2 = T(jvec+1,i,iklm1);
                    
                left = jvec < i;

                % case 3 or 4, use common formula (2.15)
                % T(j,i,n): m x 1
                jidx = loopidx & left;
                if any(jidx)
                    j = jvec(jidx);
                    m = size(j,1);
                    dt1 = t(i)-t(j);
                    dt2 = t(j+l)-t(i);
                    T1 = Tvec1(jidx);
                    T2 = Tvec2(jidx);
                    numer = dt1.*T1 + dt2.*T2;
                    denom = t(j+l)-t(j);
                    idnom = 1./denom;
                    f = numer .* idnom; % m x 1
                    T(j,i,n) = f + T(j,i+1,ikm1l);
                    
                    dti = dt(i,:); % 1 x p
                    dtj = dt(j,:); % m x p
                    dtjpl = dt(j+l,:); % m x p
                    ddt1 = bsxfun(@minus, dti, dtj);
                    ddt2 = bsxfun(@minus, dtjpl, dti);
                    dn1 = bsxfun(@times, ddt1, T1) + bsxfun(@times, ddt2, T2); % m x p
                    dn1 = reshape(dn1,[m 1 1 p]); % m x 1 x 1 x p
                    dn2 = bsxfun(@times, dt1, dT(j,i,iklm1,:)) + ...
                          bsxfun(@times, dt2, dT(j+1,i,iklm1,:)); % m x 1 x 1 x p
                    dnumer = dn1+dn2; % m x 1 x 1 x p
                    ddenom = reshape(dtjpl-dtj, [m 1 1 p]); % m x 1 x 1 x p
                    df = dnumer - bsxfun(@times, f, ddenom);
                    df = bsxfun(@times, df, idnom);
                    dT(j,i,n,:) = df + dT(j,i+1,ikm1l,:);
                end % (2.15)
                       
                % case 1, use formula (2.14)
                jidx = loopidx & ~left;
                if any(jidx)
                    j = jvec(jidx);
                    m = size(j,1);
                    dt1 = t(i+k)-t(j);
                    dt2 = t(j+l)-t(i+k);
                    T1 = Tvec1(jidx);
                    T2 = Tvec2(jidx);
                    numer = dt1.*T1 + dt2.*T2;
                    denom = t(j+l)-t(j);
                    idnom = 1./denom;
                    f = numer .* idnom; % m x 1
                    % T(j,i,n): m x 1
                    T(j,i,n) = f + T(j,i,ikm1l);
                    
                    dtipk = dt(i+k,:); % 1 x p
                    dtjpl = dt(j+l,:); % m x p
                    dtj = dt(j,:); % m x p
                    ddt1 = bsxfun(@minus, dtipk, dtj);
                    ddt2 = bsxfun(@minus, dtjpl, dtipk);
                    dn1 = bsxfun(@times, ddt1, T1) + bsxfun(@times, ddt2, T2); % m x p
                    dn1 = reshape(dn1,[m 1 1 p]); % m x 1 x 1 x p
                    dn2 = bsxfun(@times, dt1, dT(j,i,iklm1,:)) + ...
                          bsxfun(@times, dt2, dT(j+1,i,iklm1,:)); % m x 1 x 1 x p
                    dnumer = dn1+dn2; % m x 1 x 1 x p
                    ddenom = reshape(dtjpl-dtj, [m 1 1 p]); % m x 1 x 1 x p
                    df = dnumer - bsxfun(@times, f, ddenom);
                    df = bsxfun(@times, df, idnom);
                    dT(j,i,n,:) = df + dT(j,i,ikm1l,:);                   
                end % (2.14)
                
            else % k>l
                % already computed, because [K(:) L(:)] is sorted!!!
                % And T is unchanged when we swap k<->l and i<->j
                j = jvec(loopidx);
                ikl = idxkl(l+1,k+1);
                T(j,i,n) = T(i,j,ikl).';
                dT(j,i,n,:) = permute(dT(i,j,ikl,:),[2 1 3 4]);
            end % if k<=l

        end % dti~=0
    end % i-for loop
end % k/l for loop

% Extract only the part we are interested in
subsij = 1:nknots-r;
subskl = idxkl(r+1, r+1); % scalar
R = T(subsij,subsij,subskl);
dR = dT(subsij,subsij,subskl,:);
% Remove the singleton
dR = reshape(dR,[size(R) p]);

% Symmetrize
R = 0.5*(R+R.');

% Fix the bug, 08 june 2010
% dR(isnan(dR)) = 0;
dR = 0.5*(dR + permute(dR,[2 1 3]));

end % quadengine

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
function [B dB] = Btt(t, j1, j2, k, dt, leftflag)
% function [B dB] = Btt(t, j1, j2, k, dt, leftflag)
% Take the derivative with respect to t

% Linear of constant spline, nodes value stay constant to 0 or 1
% Derivatives are nil
if k<=2
    p = size(dt,2);
    B = Bernstein(t(j1), t, j2, k, [], leftflag);
    dB = zeros([size(B) p]);
else
    % B is the B-spline B_j2,k computed at t(j1)
    % Derivative with respect to the second argument knot t(j2) at direction dt
    [B dB] = BernKnotDeriv(t(j1), t, j2, k, dt, [], [], leftflag);
    m = length(j1);
    p = size(dt,2);
    
    % Derivative of B with respect to the first argument t(j1)
    [Dr td kd] = DerivB(t, k, 1, [], j2); % Dr = something x m
    B1 = Bernstein(t(j1), td, [], kd, Dr, leftflag); % m x n
    dtj1 = reshape(dt(j1,:),[m 1 p]); % m x 1 x p
    dB1 = bsxfun(@times, full(B1), dtj1); % m x n x p
    
    % Add both together
    dB = dB + dB1;
end

end % Btt

