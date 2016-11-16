function B = KramerBruckner(A, window)
%
% function B = KramerBruckner(A, window)
%
% Kramer & Bruckner filtering
%
% AUTHOR: Bruno Luong <brunoluong@yahoo.com>
% HISTORY
%   Original: 12-Jul-2009
%

if nargin<2
    window = [];
end

% cast
classA = class(A);

% Min/Max filtering
A = double(A);
[minimg maximg] = minmaxfilt(A, window, 'both', 'same');

% Kramer & Brunckner
meanimg = (minimg+maximg)/2;
B = minimg;
B(A>=meanimg) = maximg(A>=meanimg);

if ~strcmp(classA,'double')
    B = feval(classA,B);
end