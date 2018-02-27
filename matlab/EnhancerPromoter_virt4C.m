function [] = EnhancerPromoter_virt4C(fnij,fme,fgenes,res,ch,out52)

if ischar(res)
	res = str2num(res);
end;
if ischar(ch)
	ch = ['chr' ch];
else
	ch = ['chr' num2str(ch)];
end;

Log('loading genes');
[chr,from,to,gene,str]=textread(fgenes,'%s%d%d%s%[^\n]\n','delimiter','\t','headerlines',0,'bufsize',12e3);
[~,J]=sort(chr); gene=gene(J); chr=chr(J); from=from(J); to=to(J); str=str(J); clear J;
N=length(chr);
Log();

step = res;
flank=2.4e6;
pos = from; for i=1:length(pos), if length(str{i})>0 && str{i}(1)=='-', pos(i)=to(i); end; end; 
XX=-flank:step:flank;

Log('Loading matrixes');
nij = load(fnij);
hrr = load(fme);
% zero all digonal items - as we cannot find interactions with ourselves
nij(find(eye(size(nij,1))))=0;
Log();

rng('shuffle');

% output directory
if exist('../virt4C') ~= 7
    mkdir('../virt4C');
end

for i=1:N,
    if ~strcmp(chr{i},ch)
	    continue; %skip genes not in Hi-C chromosome
    end;

    fprintf('\r%d/%d %s      ',i,N,chr{i});
    [I,F,D,R,PV,FD,ZSC,X] = find_int(nij,hrr,pos(i),flank,0.05,res);

    n=length(D); mid=ceil(n/2); F(mid)=NaN; D(mid)=NaN;
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
    % legend({'Hi-C','model','Hi-C (FDR<1e-2)'},'Location','NW','FontSize',12)
    set(gcf,'PaperSize',[5 4],'PaperPosition',[0 0 5 4],'PaperPositionMode','manual');
    print(gcf,'-dpdf','-r100', sprintf('../virt4C/%s.%s.%s-%d.virt4C.pdf', gene{i}, stub, chr{i}, pos(i)));

end
fprintf('\r                   \n');
