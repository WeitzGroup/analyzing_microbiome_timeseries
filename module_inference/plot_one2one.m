function plot_one2one()

hold on;
tmpx = xlim;
tmpy = ylim;
plot(tmpx,tmpy,'k-','HandleVisibility','off');
xlim(tmpx);
ylim(tmpy);
hold off;

end