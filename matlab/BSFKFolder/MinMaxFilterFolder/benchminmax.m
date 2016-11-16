function  benchminmax()
% function  benchminmax()
%
% Benchmark runtime minmax filter among four methods: Matlab vectorized,
% for-loop, Lemire's MEX files, and vanherk by  Frederico D'Almeida
%

try
    lemire_engine(1,1);
catch
    minmaxfilter_install();
end

% Data
fprintf('Random, win=100\n');
B = randperm(100000);
win = 100;
b(B, win);

fprintf('Ordered, win=100\n');
B = -(1:100000);
win = 100;
b(B, win);

fprintf('Random, win=3\n');
B = randperm(100000);
win = 3;
b(B, win);

fprintf('Ordered, win=3\n');
B = -(1:100000);
win = 3;
b(B, win);

end % benchminmax

function b(B, win)

N = 5; % number of tests
t1 = Inf; t2 = Inf; t3 = Inf; t4=Inf;

for j = 1:N
    % Vectorized engine
    tic
    subs = ones(win*length(B)-(win-1)*win,1);
    subs(win+1:win:end) = 2-win;
    maxdata = max(reshape(B(cumsum(subs)),win,[]));
    t1=min(t1,toc);
    
    % For-loop engine
    tic
    maxdata = zeros(1,length(B)-win+1);
    m = max(B(1:win));
    for k=win+1:length(B)
        maxdata(k-win) = m;
        if B(k-win) < m
            m = max(m, B(k));
        else
            % Matt Fig's improvement
            m = B(1+k-win);
            for ii = 1+k-win+1:k
                if B(ii)>m
                    m = B(ii);
                end
            end
            % m = max(B(1+k-win:k));
        end
    end
    maxdata(end) = m;
    t2=min(t2,toc); 
    
    % Lemire's engine
    tic
    %maxdata = minmaxfilt1(B, win, 'max');
    [trash maxdata] = lemire_engine(B, win);
    t3=min(t3,toc); 
    
    % Lemire's engine
    tic
    [maxdata] = vanherk(B, win, 'max');
    t4=min(t4,toc); 
end

fprintf('\tAbsolute time [s]\n')
fprintf('\t\thankel: %f, for-loop: %f, lemire: %f, vanherk: %f\n', ...
        t1, t2, t3, t4);
tmin = min([t1 t2 t3 t4]);
fprintf('\tRelative time\n')
fprintf('\t\thankel: %f, for-loop: %f, lemire: %f, vanherk: %f\n', ...
        t1/tmin, t2/tmin, t3/tmin, t4/tmin);

end