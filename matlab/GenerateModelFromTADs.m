% GenerateModelFromTADs
%
% This beloved function creates a model from matrix and domain, and save it to disk
% Output file format will be mat file with ....
function [] = GenerateModelFromTADs(fMatrix,fDomains,res,win,fOut)
if ischar(win)
	win = str2num(win);
end;
if ischar(res)
	res = str2num(res);
end;

MIN_DIAG = 1;%0 is the middle one; %SETTING MINDIAG TO HIGHER IS DANGEROUS
MAX_DIAG = win/res; %So that 2mb in 20k window will result in such and so

%Load
a = load(fMatrix);
d = load(fDomains);

a_size = size(a,2);

%Trim diagonals
Log('Creating diag matrix');
clear_diag = CreateDiagMatrix(a_size,MIN_DIAG,MAX_DIAG);
a_diag = a .* clear_diag;
Log();

%Use public TAD data to extract intraTAD stuff
Log('Extracting TAD/Bg matrix-masks'); %STORE THIS!!!! TAKES FOREVER
[a_bgmask,a_tadmask,a_crossmask] = GetTADsMask(d,MIN_DIAG,MAX_DIAG,a_size,res);

% Depracated
%a_bgmask = a_bgmask + a_crossmask; % A convetsion so that all non-basic TAD things will be bg
%a_crossmask = a_crossmask-a_crossmask; % Zero out a_cross so it will not disturb us
% Instead, we remove the crosses because they are a mixed state

a_tad = a_tadmask .* a_diag; %No tiles (xor), just TADs
a_bg = a_bgmask .* a_diag; %Background, no cross tads
a_cross = a_crossmask .* a_diag; %Cross TAD stuff
%To calculate priors, we set a_cross to 0 so that prior(tad) = #(tad)/#(tad)+#(bg)
[pr_tad, pr_bg, pr_crs] = GetPdTBC(a_tadmask,a_bgmask,a_crossmask.*0,MIN_DIAG,MAX_DIAG);%Model
Log();

%Find gaussians
Log('Gaussianing');
a_gmm_tad = AnalyzeGMM(a_tad,MIN_DIAG,MAX_DIAG,1);
a_gmm_bg = AnalyzeGMM(a_bg,MIN_DIAG,MAX_DIAG,1);
meansTadDxn = zeros(size(a_gmm_tad));
meansBgDxn = zeros(size(a_gmm_bg));
sigmaTadDxn = zeros(size(a_gmm_tad));
sigmaBgDxn = zeros(size(a_gmm_bg));
for i = 1:numel(a_gmm_tad)
	meansTadDxn(i) = a_gmm_tad{i}.mu(1);
	meansBgDxn(i) = a_gmm_bg{i}.mu(1);
	sigmaTadDxn(i) = sqrt(a_gmm_tad{i}.Sigma(1));
	sigmaBgDxn(i) = sqrt(a_gmm_bg{i}.Sigma(1));
end
Log();

Log('Saving');
save(fOut,'meansTadDxn','meansBgDxn','sigmaTadDxn','sigmaBgDxn','pr_tad','pr_bg');
Log();

if 0 %Print Mean & Sigma for Dixon
	x_vals = (MIN_DIAG:MAX_DIAG)*res;
	figure;
	loglog(x_vals,meansTadDxn,'rx-','LineWidth',2);
	hold on;
	loglog(x_vals,meansBgDxn,'bx-','LineWidth',2);
	legend('Tad','Bound');

	loglog(x_vals,meansTadDxn+sigmaTadDxn,'r--');
	loglog(x_vals,meansTadDxn-sigmaTadDxn,'r--');
	loglog(x_vals,meansBgDxn+sigmaBgDxn,'b--');
	loglog(x_vals,meansBgDxn-sigmaBgDxn,'b--');
	title('Mean vs. Distance');
	xlabel('Distance (bp)');
	ylabel('Mean');
	fprintf('Done\r\n');
end

%Plot PdT, pdB, pdC
if 0
	x_vals = (MIN_DIAG:MAX_DIAG)*res;
	figure; hold on;
	plot(x_vals,pr_tad,'LineWidth',2);
	plot(x_vals,pr_bg,'LineWidth',2);
	legend('TAD','Bg');
	xlabel('Distance (bp)');
	ylabel('Prob');
	title('Pr(TAD), Pr(Bg) vs. Distance');
end

end
