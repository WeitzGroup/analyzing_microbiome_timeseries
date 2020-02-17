function plot_timeseries(t,X,host_or_virus,fs)

if strcmp(host_or_virus,'host')
    titlestr = 'H (host samples)';
    ystr = 'hosts';
    cmapstr = 'parula';
else
    titlestr = 'W (transformed virus differences)';
    ystr = 'viruses';
    cmapstr = 'hot';
end

%t = t(1:end-1);
imagesc(X);
title(titlestr);
set(gca,'YTick',[]);
ylabel(ystr);
tmpID = round(linspace(1,length(t),5));
set(gca,'XTick',tmpID-0.5,'XTickLabel',string(t(tmpID)));
xlabel('time (hrs)');
colormap(gca,cmapstr);
colorbar();
set(gca,'FontSize',fs);

end