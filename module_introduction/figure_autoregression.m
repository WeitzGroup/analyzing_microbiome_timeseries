
if ~exist('X','var')
    fprintf('running autoregression tutorial...\n');
    tutorial_autoregression
end

fs = 8;
fig = figure();
nrow = 3;
ncol = 2;

subplot(nrow,ncol,1);
plot(X);
xlabel('time');
ylabel('X');
title('Independent random walks');
set(gca,'FontSize',fs);

subplot(nrow,ncol,2);
plot_correlation(rho_X,pval_X);
title(sprintf('Correlations among\nindependent random walks'));
set(gca,'FontSize',fs);

ax = subplot(nrow,ncol,3);
plot(Xresidual);
tmpy = ylim;
ylim([-1 +1]*max(abs(tmpy)));
xlabel('time');
ylabel('X_{residual}');
title(sprintf('Residuals of\nindependent random walks'));
set(gca,'FontSize',fs);
lg = legend(string(1:N),'FontSize',6,'Location','eastoutside','NumColumns',3);
lg.Title.String = 'timeseries #';
lg.Position = [ax.Position(1) ax.Position(2)-0.15 ax.Position(3) 0.05];

subplot(nrow,ncol,4);
plot_correlation(rho_Xresidual,pval_Xresidual);
title(sprintf('Correlations among\nresidual timeseries'));
set(gca,'FontSize',fs);

fig.Units = 'inches';
fig.Position(3:4) = [6 6];
print('figure_autoregression','-djpeg','-r300');