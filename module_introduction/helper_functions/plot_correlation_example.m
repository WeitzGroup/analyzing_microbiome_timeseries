function plot_correlation_example(X1,X2,original_or_residual)
% Plot timeseries X1 vs X2

subplot(1,2,1);
plot(X1);
hold on;
plot(X2);
hold off;
xlabel('time');
ylabel(sprintf('X %s',original_or_residual));
legend({'X1','X2'},'Location','SouthEast');
set(gca,'FontSize',14);

subplot(1,2,2);
plot(X1,X2,'o');
xlabel(sprintf('X1 %s',original_or_residual));
ylabel(sprintf('X2 %s',original_or_residual));
text(0.95,0.1,sprintf('r=%.2f',corr(X1,X2)),...
    'Units','normalized','HorizontalAlignment','right',...
    'FontSize',16,'FontWeight','bold','Color','r');
title(sprintf('Example of a "highly correlated" pair (%s timeseries)',original_or_residual),...
    'HorizontalAlignment','right');
set(gca,'FontSize',14);

tmppos = get(gcf,'Position');
set(gcf,'Position',tmppos.*[1 1 1.5 .75]);

end