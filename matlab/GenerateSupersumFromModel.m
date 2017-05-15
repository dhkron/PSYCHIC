% GenerateProbmapsFromModel
%
% model+matrix to probTAD,probBG and LPR 
function [] = GenerateSupersumFromModel(fMatrix,fModel,res,win,fOutSupersum,fOutT,fOutB,fOutLLR)
if ischar(win)
	win = str2num(win);
end;
if ischar(res)
	res = str2num(res);
end;

MIN_DIAG = 1;%0 is the middle one; %SETTING MINDIAG TO HIGHER IS DANGEROUS
MAX_DIAG = win/res; %So that 2mb in 20k window will result in such and so

%Load
Log('Loading')
a = load(fMatrix)/1;
m = load(fModel);
Log()

a_size = size(a,2);

Log('Generating Matrixes')
%True for not performing per-cell normalization
[a_llr,a_pdt,a_pdb] = PosteriorHeatmap(a,MIN_DIAG,MAX_DIAG,m.pr_bg,m.meansBgDxn,m.sigmaBgDxn,m.pr_tad,m.meansTadDxn,m.sigmaTadDxn,false);
%[a_llr,a_pdt,a_pdb] = LiklihoodHeatmap(a,MIN_DIAG,MAX_DIAG,m.meansBgDxn,m.sigmaBgDxn,m.meansTadDxn,m.sigmaTadDxn);
a_pdt_log = log(a_pdt);
a_pdt_log(isinf(a_pdt_log)) = 0;
a_pdb_log = log(a_pdb);
a_pdb_log(isinf(a_pdb_log)) = 0;
Log()

Log('Supersumming')
%Sum of logs accross pyramids
%supersum_tad = PyramidSky(a_pdt_log,0,MAX_DIAG,MIN_DIAG);
%supersum_bg = PyramidSky(a_pdb_log,0,MAX_DIAG,MIN_DIAG);
%supersum = supersum_tad-supersum_bg;
supersum = PyramidSky(a_llr,-a_llr,MAX_DIAG,MIN_DIAG);
Log()

Log('Saving')
dlmwrite(fOutT,a_pdt);
dlmwrite(fOutB,a_pdb);
dlmwrite(fOutSupersum,supersum);
dlmwrite(fOutLLR,a_llr);
Log()

end
