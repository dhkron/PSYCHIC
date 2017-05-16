function [] = virt4C(name,chr,pos,nij,hrr,flank,res)
% % for example, to create a virt4C plot for Foxg1 (1Mb flanks, 40Kb Hi-C resolution):
% name = 'Foxg1'; chr='chr12'; pos=50480000; flank=1e6; res=40e3;
% nij=load('../examples/output/mCO.chr12.matrix.txt'); nij(find(eye(size(nij,1))))=0;
% hrr=load('../examples/output/mCO.chr12.model.estimated.matrix.txt');

% name = 'Shh'; chr='chr5', pos=28800000; flank=1e6; res=40e3;

[inter,F,D,R,PV,FD,V,X] = find_int(nij,hrr,pos,1e6,1e-3,res);
n=length(D); mid=ceil(n/2); F(mid)=NaN; D(mid+[-1:1])=NaN;

clrs=jet(10);
bar(X,D,.5,'FaceColor',clrs(4,:),'EdgeColor','none');
hold on;
plot(X,F,'k-','LineWidth',4);
bar(X,(FD<=1e-2).*D,.5,'FaceColor',clrs(10,:),'EdgeColor','none');
hold off;

axis tight;
XX=40e3*(X-X(1))/1e6; XX=XX-XX(end)/2;
II=find(ceil(XX)==XX);

tmpstr = sprintf('%dMb,',XX(II));
xt = strsplit(tmpstr(1:end-1),','); xt{ceil(length(xt)/2)}=name;

set(gca,'YLim',[.1 20],'XTick',X(II),'XTickLabel',xt,'FontSize',16);
title('Virtual 4C'); ylabel('Intensity');
legend({'Hi-C','Background model','Hi-C (FDR<1e-2)'},'Location','NW')
set(gcf,'PaperSize',[5 4],'PaperPosition',[0 0 5 4],'PaperPositionMode','manual');
print(gcf,'-dpdf','-r100', sprintf('../examples/%s.%s-%d.virt4C.pdf', name, chr, pos));

