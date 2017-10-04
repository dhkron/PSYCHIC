function [] = EnhancerPromoter_virt4C(fnij,fme,fgenes,res,ch,out4,out3,out2,outRand,out52)
%% modify EnhancerPromoter.m to generate virtual 4C plots
%% for a specific chromosome.
%% Only generates for one gene (line 52)

if ischar(res)
	res = str2num(res);
end;
if ischar(ch)
	ch = ['chr' ch];
else
	ch = ['chr' num2str(ch)];
end;

Log('loading genes');
[chr,from,to,gene,~,str]=textread(fgenes,'%s%d%d%s%d%c%*[^\n]\n','delimiter','\t','headerlines',0,'bufsize',12e3);
[~,J]=sort(chr); gene=gene(J); chr=chr(J); from=from(J); to=to(J); str=str(J); clear J;
N=length(chr);
Log();

step = res;
flank=1.0e6;
pos = from; for i=1:length(pos), if str(i)=='-', pos(i)=to(i); end; end;
XX=-flank:step:flank;

Log('Loading matrixes');
nij = load(fnij);
hrr = load(fme);
% zero all digonal items - as we cannot find interactions with ourselves
nij(find(eye(size(nij,1))))=0;
Log();

virt4C = 1;
if ~virt4C
    fid1 = fopen(out4,'w');
    fid2 = fopen(out3,'w');
    fid3 = fopen(out2,'w');
    fid4 = fopen(outRand,'w');
    fid5 = fopen(out52,'w');
end
rng('shuffle');

for i=1:N,
    if ~strcmp(chr{i},ch)
	    continue; %skip genes not in Hi-C chromosome
    end;

    fprintf('\r%d/%d %s      ',i,N,chr{i});
    [I,F,D,R,PV,FD,ZSC,X] = find_int(nij,hrr,pos(i),flank,0.05,res);

    if virt4C,
	if ~strcmpi(gene{i},'SHH'), continue; end;

	n=length(D); mid=ceil(n/2); F(mid)=NaN; D(mid+[-1:1])=NaN;

	clrs=jet(10);
	bar(X,D,.5,'FaceColor',clrs(4,:),'EdgeColor','none');
	hold on;
	plot(X,F,'k-','LineWidth',4);
	bar(X,(FD<=1e-2).*D,.5,'FaceColor',clrs(10,:),'EdgeColor','none');
	hold off;

	axis tight; ax=axis;
	XX = (res*(X-X(1)) - flank ) / 1e6;
	II=find(ceil(XX)==XX);

	tmpstr = sprintf('%dMb,',XX(II));
	xt = strsplit(tmpstr(1:end-1),','); xt{ceil(length(xt)/2)}=gene{i};

	ii=max(strfind(['#',out52],'/')); jj=min(strfind(out52(ii:end),'.')); stub=out52(ii:ii+jj-2);

	set(gca,'YLim',[.1 ax(4)],'XTick',X(II),'XTickLabel',xt,'FontSize',16);
	title(sprintf('%s (%s)',gene{i},stub)); ylabel('Intensity');
	% legend({'Hi-C','Background model','Hi-C (FDR<1e-2)'},'Location','NW')
	set(gcf,'PaperSize',[5 4],'PaperPosition',[0 0 5 4],'PaperPositionMode','manual');
	print(gcf,'-dpdf','-r100', sprintf('../virt4C/%s.%s.%s-%d.virt4C.pdf', gene{i}, stub, chr{i}, pos(i)));

	break;
    end
	break;

    G(i).name = gene{i}; G(i).chr = chr{i}; G(i).pos = pos(i); G(i).hic = D; G(i).resi = R; G(i).pval = PV; G(i).fit = F; G(i).fdr = FD; 

    ppos=pos(i); gene_ind = ceil(ppos/step);
    
    I=find(FD<=1e-4);
    for j=[I]
        enh_end = step*gene_ind + XX(j);
        % fprintf(fid1,'%s\t%d\t%d\t%s:%.0fKb:%.2g\n', chr{i}, enh_end-step, enh_end, G(i).name, XX(j)/1e3, FD(j));
        fprintf(fid1,'%s\t%d\t%d\t%s:%.0fKb:%.2g:%.2g:%.2g:%.2g\n', chr{i}, enh_end-step, enh_end, G(i).name, XX(j)/1e3, FD(j), PV(j), F(j), D(j));
    end

    I=find(FD<=1e-3);
    for j=[I]
        enh_end = step*gene_ind + XX(j);
        % fprintf(fid2,'%s\t%d\t%d\t%s:%.0fKb:%.2g\n', chr{i}, enh_end-step, enh_end, G(i).name, XX(j)/1e3, FD(j));
	fprintf(fid2,'%s\t%d\t%d\t%s:%.0fKb:%.2g:%.2g:%.2g:%.2g\n', chr{i}, enh_end-step, enh_end, G(i).name, XX(j)/1e3, FD(j), PV(j), F(j), D(j));
    end
    
    I=find(FD<=1e-2);
    for j=[I]
        enh_end = step*gene_ind + XX(j);
        % fprintf(fid3,'%s\t%d\t%d\t%s:%.0fKb:%.2g\n', chr{i}, enh_end-step, enh_end, G(i).name, XX(j)/1e3, FD(j));
	fprintf(fid3,'%s\t%d\t%d\t%s:%.0fKb:%.2g:%.2g:%.2g:%.2g\n', chr{i}, enh_end-step, enh_end, G(i).name, XX(j)/1e3, FD(j), PV(j), F(j), D(j));
    end
    
    I=find(rand(1,length(PV))<1e-2);
    PVrnd = rand(1,length(PV)); FDrnd = PVrnd;
    % [~,FDrnd] = mafdr(PVrnd,'Lambda',0.05);
    I=find(PVrnd<=1e-2);
    for j=[I]
        enh_end = step*gene_ind + XX(j);
        % fprintf(fid4,'%s\t%d\t%d\t%s:%.0fKb:%.2g\n', chr{i}, enh_end-step, enh_end, G(i).name, XX(j)/1e3, FD(j));
	fprintf(fid4,'%s\t%d\t%d\t%s:%.0fKb:%.2g:%.2g:%.2g:%.2g\n', chr{i}, enh_end-step, enh_end, G(i).name, XX(j)/1e3, FDrnd(j), PVrnd(j), F(j), D(j));
    end

    I=find(FD<=5e-2);
    for j=[I]
        enh_end = step*gene_ind + XX(j);
        % fprintf(fid5,'%s\t%d\t%d\t%s:%.0fKb:%.2g\n', chr{i}, enh_end-step, enh_end, G(i).name, XX(j)/1e3, FD(j));
	fprintf(fid5,'%s\t%d\t%d\t%s:%.0fKb:%.2g:%.2g:%.2g:%.2g\n', chr{i}, enh_end-step, enh_end, G(i).name, XX(j)/1e3, FD(j), PV(j), F(j), D(j));
    end

end
fprintf('\r                   \n');
fclose(fid1); fclose(fid2); fclose(fid3); fclose(fid4); fclose(fid5);
