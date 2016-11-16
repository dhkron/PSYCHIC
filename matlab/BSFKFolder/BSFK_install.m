function BSFK_install
% function BSFK_install
% Install BSFK package
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% History: Original 17-Nov-2009
%          19-Nov-2009: save the path
%          04-Dec-2009: minmaxfilt package

disp('Installing BSFK');
path = fileparts(mfilename('fullpath'));

newpath = path();
disp(['addpath(''' newpath ''')']);
addpath(newpath);

if exist('MultiSolver','file') < 2
    newpath = [path filesep 'MultiSolverFolder'];
    disp(['addpath(''' newpath ''')']);
    addpath(newpath);
end

if exist('pseudoinverse','file') < 2
    newpath = [path filesep 'PseudoInverseFolder'];
    disp(['addpath(''' newpath ''')']);
    addpath(newpath);
end

if exist('minmaxfilt','file') < 2
    newpath = [path filesep 'MinMaxFilterFolder'];
    disp(['addpath(''' newpath ''')']);
    addpath(newpath);
    currpath = cd(newpath);
    minmaxfilter_install();
    cd(currpath);
end

if exist('quadprog','file') < 2
    disp('---------------------------------------------------------------')
    disp('Additional steps: download and install either QP engine at')
    fprintf('\t http://sigpromu.org/quadprog/index.html or\n');
    fprintf('\t http://www.mat.univie.ac.at/~neum/software/minq/\n');
    disp('---------------------------------------------------------------')
end
savepath();
    