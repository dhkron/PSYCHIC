function minmaxfilter_install
% function minmaxfilter_install
% Installation by building the C-mex files for minmax filter package
%
% Author Bruno Luong <brunoluong@yahoo.com>
% Last update: 20-Sep-2009

arch=computer('arch');
mexopts = {'-O' '-v' ['-' arch]};
% 64-bit platform
if ~isempty(strfind(computer(),'64'))
    mexopts(end+1) = {'-largeArrayDims'};
end

% invoke MEX compilation tool
mex(mexopts{:},'lemire_nd_minengine.c');
mex(mexopts{:},'lemire_nd_maxengine.c');

% Those mex files are no longer used, we keep for compatibility
mex(mexopts{:},'lemire_nd_engine.c');
mex(mexopts{:},'lemire_engine.c');