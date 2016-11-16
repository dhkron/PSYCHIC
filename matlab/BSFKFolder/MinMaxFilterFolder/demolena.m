% demo KramerBruckner filtering

try
    lenaorg=imread('lena_std.tif');
catch
    site = 'http://www-2.cs.cmu.edu/~chuck/lennapg/lena_std.tif';
    fprintf('Download first lena_std.tif from <%s>\n', site);
    fprintf('then rerun the script\n');
    return
end

try
    lemire_engine(1,1);
catch
    minmaxfilter_install();
end

lenaorg=sum(double(lenaorg),3);

lena = lenaorg;
for k=1:10
    lena = KramerBruckner(lena, [15 15]);
end

figure(1);
ax=subplot(1,2,1,'Parent',1);
imagesc(lenaorg,'Parent',ax);
axis(ax,'equal')
set(ax,'Visible','off');
colormap(ax,'gray')
ax=subplot(1,2,2,'Parent',1);
set(ax,'Visible','off');
imagesc(lena,'Parent',ax);
axis(ax,'equal')
set(ax,'Visible','off');
colormap(ax,'gray')