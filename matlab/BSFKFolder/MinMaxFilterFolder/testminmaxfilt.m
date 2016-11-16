% Author Bruno Luong <brunoluong@yahoo.com>
% Last update: 20-Sep-2009
% Script to tets min/max filtering

clear

try
    lemire_engine(1,1);
catch
    minmaxfilter_install();
end

fprintf('Test started... ');

%%%%%%%%%
% 1D filtering
n = 100;
w = 30;
for ntest=1:100
    A=uint8(ceil(rand(1,n)*100));
    [minA maxA] = minmaxfilt1(A, w, 'both', 'full');
    minAc = zeros(size(minA));
    maxAc = zeros(size(minA));
    for k=1:n+(w-1)
        Awin = A(max(k-w+1,1):min(k,end));
        minAc(k) = min(Awin);
        maxAc(k) = max(Awin);
    end
    if ~isequal(minAc,minA) || ~isequal(maxAc,maxA)
        fprintf('Something is wrong in 1D filter\n');
        keyboard;
    end        
end

clear

%%%%%%%%
% 2D filtering
A=rand(1024);
win=[3 5];
tic
[minA maxA minidx maxidx] = minmaxfilt(A,win,'both','full');
toc

if ~isequal(minA,A(minidx)) || ~isequal(maxA,A(maxidx))
    fprintf('Something is wrong in ND indexing\n');
    keyboard;
end
    

for i=1:size(A,1)+win(1)-1
    idxi = max(1-win(1)+i,1):min(i,size(A,1));
    for j=1:size(A,2)+win(2)-1
        idxj = max(1-win(2)+j,1):min(j,size(A,2));
        Aij = A(idxi,idxj);
        a = min(Aij(:));
        b = max(Aij(:));
        if a~=minA(i,j) || b~=maxA(i,j)
            fprintf('Something is wrong in ND filter\n');
            keyboard
        end
    end
end

fprintf('Congratulation, the minmax filter successfully tested\n');