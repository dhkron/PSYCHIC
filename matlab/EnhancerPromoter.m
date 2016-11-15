function [] = EnhancerPromoter(fnij,fme,fgenes,res,ch,out4,out3,out2,outRand,out52)
% ch is passed for filtering the gene file

if ischar(res)
	res = str2num(res);
end;
if ischar(ch)
	ch = ['chr' ch];
else
	ch = ['chr' num2str(ch)];
end;

Log('loading genes');
[chr,from,to,gene,~,str]=textread(fgenes,'%s%d%d%s%d%c%*[^\n]\n','delimiter','\t','headerlines',0,'bufsize',12e3); %fgenes was mm9.genes2.bed
[~,J]=sort(chr); gene=gene(J); chr=chr(J); from=from(J); to=to(J); str=str(J); clear J;
N=length(chr);
Log();

step = res; %was 40000
flank=1e6;
pos = from; for i=1:length(pos), if str(i)=='-', pos(i)=to(i); end; end;
XX=-flank:step:flank;

Log('Loading matrixes');
nij = load(fnij);
hrr = load(fme);
nij(find(eye(size(nij,1))))=0; %Moved from find_int.m.org
Log();

fid1 = fopen(out4,'w'); %out4 was 'mCO_1e-4.bed'
fid2 = fopen(out3,'w'); %out3 was 'mCO_1e-3.bed'
fid3 = fopen(out2,'w'); %out2 was 'mCO_1e-2.bed'
fid4 = fopen(outRand,'w'); %outRand was 'mCO_rand.bed'
fid5 = fopen(out52,'w');

rng('shuffle');

for i=1:N,
    if ~strcmp(chr{i},ch)
	    continue; %skip genes not in chromosome
    end;
    fprintf('\r%d/%d %s      ',i,N,chr{i});
    [I,FDR,D,R,PV,FD,ZSC] = find_int(nij,hrr,pos(i),flank,0.05,res); %flank was 1e6
    G(i).name = gene{i};
    G(i).chr = chr{i};
    G(i).pos = pos(i);
    G(i).hic = D;
    G(i).resi = R;
    G(i).pval = PV;
    G(i).fdr = FD;

    % ppos=pos(i)-mod(pos(i),step)+step/2 % should be correct
    % ppos=pos(i)-step/2; % what I did before (off)
    ppos=pos(i); % gives the best results, a bit shifted to the left
    I=find(FD<=1e-4);
    for j=[I]
	fprintf(fid1,'%s\t%d\t%d\t%s:%.0fKb:%.2g\n', chr{i}, ppos+XX(j)-step, ppos+XX(j), G(i).name, XX(j)/1e3, FD(j));
    end

    I=find(FD<=1e-3);
    for j=[I]
        fprintf(fid2,'%s\t%d\t%d\t%s:%.0fKb:%.2g\n', chr{i}, ppos+XX(j)-step, ppos+XX(j), G(i).name, XX(j)/1e3, FD(j));
    end
    
    I=find(FD<=1e-2);
    for j=[I]
        fprintf(fid3,'%s\t%d\t%d\t%s:%.0fKb:%.2g\n', chr{i}, ppos+XX(j)-step, ppos+XX(j), G(i).name, XX(j)/1e3, FD(j));
    end
    
    I=find(rand(1,length(PV))<1e-2);
    for j=[I]
        fprintf(fid4,'%s\t%d\t%d\t%s:%.0fKb:%.2g\n', chr{i}, ppos+XX(j)-step, ppos+XX(j), G(i).name, XX(j)/1e3, FD(j));
    end

    I=find(FD<=5e-2);
    for j=[I]
        fprintf(fid5,'%s\t%d\t%d\t%s:%.0fKb:%.2g\n', chr{i}, ppos+XX(j)-step, ppos+XX(j), G(i).name, XX(j)/1e3, FD(j));
    end

end
fprintf('\r                   \n');
fclose(fid1); fclose(fid2); fclose(fid3); fclose(fid4);

% save('mCO.PromoterEnhancer.mat','G');
