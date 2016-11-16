function [x,y,ier,lbd] = qpmin(c,G,A,b,eq,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% qpmin.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x,y,ier,lbd] = qpmin(c,G,A,b,eq,options)
% minimizes a definite quadratic form subject to linear constraints
%    min    fct = 0.5 x'*G*x + c'*x
%    s.t.   A*x >= b
% where G is a definite symmetric (n x n) matrix
%
% if A is sparse, it is assumed that the ordering is such that
% sparse Cholesky factorization of G and AG^(-1)A^T are feasible
%
% When the METHOD is QPAS, LBD the structure of Lagrange multipliers
% is returned:
%    lbd.equality: Lagrange multipliers for equality constraints 
%    lbd.inequality: Lagrange multipliers for general inequality constraints 
%    lbd.lowerbound: Lagrange multipliers for lower bound constraints 
%    lbd.upperbound: Lagrange multipliers for upper bound constraints
%

NOTFEASIBLE = 1;

if nargin<6 || isempty(options)
    options = struct();
end

engine = getoption(options,'engine','qpas');
prt = getoption(options, 'prt', 1);

if isempty(A)
    A = zeros(0,length(G));
    b = zeros(0,1);
end

% Detect singular modes of G, use QR rather SVD because it can operates
% on sparse matrix
[m n] = size(A);

% Initialize the eq array
if isempty(eq)
    % By default, all are inequalities
    eq = false(m,1);
elseif ~islogical(eq)
    eq = logical(eq);
end

switch lower(engine)
    case {'qpas'}
        
        [x ier lbd] = qpas(full(G),full(c),...
                          -full(A(~eq,:)),-full(b(~eq,:)),...
                          full(A(eq,:)),full(b(eq,:)),[],[],0);
                      
        y = lbd.inequality;
        
    case {'blform', 'qpaspcx'}
        
        % G(:,p) = Q*R
        % Q'*Q = I
        % R triangular
        [Q R p] = qr(full(G),0); %#ok
        clear Q
        
        % reverse the permutation
        ip = zeros(size(p),'double');
        ip(p) = 1:length(p);
        
        wrn = warning('off','MATLAB:normest:notconverge');
        deftol = n * normest(G) * eps(class(G));
        warning(wrn.state,'MATLAB:normest:notconverge');
        % Default cutting rank
        tol = getoption(options,'ranktol',deftol);
        % estimate the rank
        k = find(abs(diag(R))<=tol,1,'first');
        
        xx = getoption(options,'xguess',zeros(n,1));
        
        
        % k is the rank of G
        if k==0
            error('G matrix cannot be all zero');
        else
            if ~isempty(k)
                rankG = k-1;
                % SVD, because G is symmetric, U is equal to V
                [U,S,V] = svd(full(G),0); %#ok
                clear V
                
                % G = U*S*U'
                s = diag(S);
                clear S;
                idx1 = 1:rankG;
                idx2 = rankG+1:n;
                s(idx2) = 0;
                c = (c.'*U).'; % c = U'*c;
                d = c(idx1);
                e = c(idx2);
                H = sparse(idx1,idx1,s(idx1),rankG,rankG); %G(idx1,idx1);
                P = A*U(:,idx1);
                Q = A*U(:,idx2);
                Qt = Q.';
                lambda0 = pseudoinverse(Qt)*e; % same dimension as the constraints
                f = d-P.'*lambda0;
                xx = (xx.'*U(:,idx1)).';
                lobnd = -lambda0;
                R = null(full(Qt));
            else % full rank
                rankG = n;
                H = G;
                P = A;
                Q = zeros(size(A,1),0);
                f = c;
                lobnd = zeros(m,1);
                %xx = xx;
            end
        end
        %
        
        % Number of constraints
        m = size(P,1);
      
        % Matrix and rhs for dual problem
        Gdual = P*(H\(P.'));
        cdual = -(P*(H\f) + b);
        
        % there is no lower bound for equality constraint
        lobnd(eq) = -inf;
        
        if rankG == n
            L = [];
            k = [];
            l = lobnd;
        else
            Gdual = R.'*Gdual*R;
            cdual = R.'*cdual;
            L = -R;
            k = -lobnd;
            l = [];
            u = [];
        end
        
        % Gdual is full rank, we can call the quadratic minimization engine
        [mu,ier,lm] = qpas(Gdual,cdual,L,k,[],[],l,u,1);
        
        if rankG == n
            xi = mu;
        else
            xi = R*mu;
        end
        
        % Recover the primal
        y = H\(P.'*xi-f);
        if ier==99, return; end;
        
        if rankG==n
            x = y;
            return
        end
        % Solve linear programming
        if ~isempty(Q)
            cost = e;
            Ac = -Q;
            bc = P*y - b;
            Aineq = sparse(Ac(~eq,:));
            bineq = bc(~eq,:);
            Aeq = sparse(Ac(eq,:));
            beq = bc(eq,:);
            % Minimizing the dual problems
            [z f LPflag] = PCx(cost, Aineq, bineq, Aeq, beq, [], []);
            
            if length(z)~= length(cost)
                z = zeros(size(cost));
            end
            
            x = [y; z];
            
            % Back to the original basis
            x = U*x;
        end
        
    otherwise
        
        error('qpmin: unknown methods')
        
end

end % qpmin
