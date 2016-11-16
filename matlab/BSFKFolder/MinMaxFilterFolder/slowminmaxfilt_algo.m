function [minval maxval] = slowminmaxfilt_algo(a, window)
%
% function [minval maxval] = slowminmaxfilt_algo(a, window)
%
% This function uses to illustrate how Lemire's streaming min/max
% works. It is NOT recommended to use in your application.
% This algorithm is used by this package
%
% Reference: Lemire's "STREAMING MAXIMUM-MINIMUM FILTER USING NO MORE THAN
% THREE COMPARISONS PER ELEMENT" Nordic Journal of Computing, Volume 13,
% Number 4, pages 328-339, 2006.
%
% AUTHOR: Bruno Luong <brunoluong@yahoo.com>
% HISTORY
%   Original: 12-Jul-2009

n=length(a);

minval = nan(1,n-window+1);
maxval = nan(1,n-window+1);

L = initwedge(window+1);
U = initwedge(window+1);

for i = 2:n
    if i > window
        if ~wedgeisempty(U)
            maxval(i-window) = a(getfirst(U));
        else
            maxval(i-window) = a(i-1);
        end
        if ~wedgeisempty(L)
            minval(i-window) = a(getfirst(L));
        else
            minval(i-window) = a(i-1);
        end
    end % i>window
    
    if a(i) > a(i-1)
        L = pushback(L, i-1);
        if i==window+getfirst(L), L = popfront(L); end
        while ~wedgeisempty(U)
            if a(i)<=a(getlast(U))
                if i==window+getfirst(U), U = popfront(U); end
                break
            end
            U = popback(U);
        end % while loop       
    else
        U = pushback(U, i-1);
        if i==window+getfirst(U), U = popfront(U); end
        while ~wedgeisempty(L)
            if a(i)>=a(getlast(L))
                if i==window+getfirst(L), L = popfront(L); end
                break
            end
            L = popback(L);
        end % while loop
    end % if
    
    fprintf('---- i=%d \n', i);
    dispw(L)
    dispw(U)

end % for-loop

i=n+1;
if ~wedgeisempty(U)
    maxval(i-window) = a(getfirst(U));
else
    maxval(i-window) = a(i-1);
end
if ~wedgeisempty(L)
    minval(i-window) = a(getfirst(L));
else
    minval(i-window) = a(i-1);
end

% L.mxn
% U.mxn

end

function X = initwedge(sz)
X = struct('buffer', zeros(1, sz), ...
           'sz', sz, ...
           'n', 0, ...
           'first', 1, ...
           'last', 0, ...
           'mxn', 0);
end

function b = wedgeisempty(X)
b = X.n <= 0;
end

function X = pushback(X, v)
X.last = mod(X.last,X.sz) + 1;
X.buffer(X.last) = v;
X.n = X.n+1;
X.mxn = max(X.mxn,X.n); % track the maximul of buffer size
end

function X = popback(X)
X.n = X.n-1;
X.last = mod(X.last-2,X.sz) + 1;
end

function X = popfront(X)
X.n = X.n-1;
X.first = mod(X.first,X.sz) + 1;
end

function v = getfirst(X)
v = X.buffer(X.first);
end

function v = getlast(X)
v = X.buffer(X.last);
end

function dispw(X)
A = X.buffer(mod(X.first+(-1:X.n-2), X.sz) + 1);
disp(A)
end