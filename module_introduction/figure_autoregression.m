
% generate figure in manuscript
% for panels A,B,D,E:
%   timeseries are newly generated every time the tutorial is run,
%   so individual timeseries and correlations will look different
% for panels C,F:
%   ensemble results are generated using the script randomwalk_ensemble
%   those data are saved, and so those panels will look the same
%   the ensemble can be rerun if desired (results are robust)

if ~exist('X','var')
    fprintf('running autoregression tutorial...\n');
    tutorial_autoregression;
    close all;
end

fs = 7;
fig = figure();
nrow = 3;
ncol = 3;

subplot(nrow,ncol,1);
plot(X);
xlabel('time');
ylabel('X');
title('Independent random walks');
set(gca,'FontSize',fs);
label_ax('A)');

subplot(nrow,ncol,2);
plot_correlation(rho_X,pval_X);
title(sprintf('Correlations among\nindependent random walks'));
set(gca,'FontSize',fs);
label_ax('B)');

subplot(nrow,ncol,3);
plot_distribution('randomwalk_ensemble/randomwalk_ensemble_originals');
title('Distribution of correlation values');
set(gca,'FontSize',fs);
label_ax('C)');

ax = subplot(nrow,ncol,4);
plot(Xresidual);
tmpy = ylim;
ylim([-1 +1]*max(abs(tmpy)));
xlabel('time');
ylabel('X_{residual}');
title(sprintf('Residuals of\nindependent random walks'));
set(gca,'FontSize',fs);
label_ax('D)');
clear tmpy;

lg = legend(string(1:N),'FontSize',6,'Location','eastoutside','NumColumns',3);
lg.Title.String = 'timeseries #';
lg.Position = [ax.Position(1) ax.Position(2)-0.15 ax.Position(3) 0.05];
clear lg ax;

subplot(nrow,ncol,5);
plot_correlation(rho_Xresidual,pval_Xresidual);
title(sprintf('Correlations among\nresidual timeseries'));
set(gca,'FontSize',fs);
label_ax('E)');

subplot(nrow,ncol,6);
plot_distribution('randomwalk_ensemble/randomwalk_ensemble_residuals');
title('Distribution of correlation values');
set(gca,'FontSize',fs);
label_ax('F)');

fig.Units = 'inches';
fig.Position(3:4) = [9 6.5];
print('figure_autoregression','-djpeg','-r300');
clear fs ncol nrow fig;