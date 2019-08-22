function plot_distribution(file_str)

load(file_str,'rho');

% histogram
histogram(rho,50,'BinLimits',[-1 +1],'Normalization','pdf');

% quantiles
Q = quantile(rho,[0.05 0.95]);

% plot
hold on;
tmpy = ylim*1.1;
plot([Q(1) Q(1)],tmpy,'r-','LineWidth',1.5);
plot([Q(2) Q(2)],tmpy,'r-','LineWidth',1.5);
text(Q(2)+0.05,tmpy(2)*0.8,sprintf('r_p=%.2f',Q(2)),'Color','red','FontSize',7);
hold off;
ylim(tmpy);
xlabel('correlation');
ylabel('PDF');

end