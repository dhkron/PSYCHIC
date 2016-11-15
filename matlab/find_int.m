function [inter,F,D,R,PV,FD,V] = find_int(nij,hrr,pos,flank,thr,res)
% inter - list of positions for interactions (relative distances)
% F - Fitted value for each position (in above list)
% D - Hi-C coverage for each window (currently 40K) in interaction with the promoter
% R - Residuals
% PV - p-values for each window in interaction
% FD - FDR for each window in interaction
% V - Z-scores for each window

inter = []; F = []; D = []; R = []; PV = []; FD = []; V=[];

% get relative positions
from=pos-flank; to=pos+flank;
step=res; %was 40000
X=max(1,floor(from/step)):min(floor(to/step),size(nij,1));

% extract data in triangle
D1=nij(X,X); F1=hrr(X,X); MM = D1 - F1; 
k2=length(X); k=ceil(k2/2);

% stop if too small (e.g. ends of chromosmoes)
if k<10, return; end;

% calculate z-scores for entire triangle
PP=NaN*MM; I=~isnan(MM); PP1=zscore(MM(I)); PP(I)=PP1; clear PP1 I;
% PP=zscore(MM(:)); PP=reshape(PP,k2,k2);
% get all windows that interact with anchor (promoter)
V=[PP(1:k,k)' PP(k,(k+1):end)];
% calculate p-values using standard Normal distribution
PV=1-normcdf(V,0,1);
% get the actual Hi-C data for these positions
D=[D1(1:k,k)' D1(k,(k+1):end)];
% and their fit
F=[F1(1:k,k)' F1(k,(k+1):end)];
% and residuals
R=D-F;

% calculate FDR value from p-values
FD = mafdr(PV,'Lambda',thr);
% find significant interactions
I=find(FD<thr);

% and extract the FDR for them
F=FD(I);
% and their relative distance from the promoter
inter = step*(I-k);
