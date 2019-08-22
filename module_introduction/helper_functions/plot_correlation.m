function plot_correlation(rho, pval, flag_labels)
% Plot correlation matrix and mark significant values.

% if flag_labels is passed, turns off the x and y tick labels

sig = 0.05; % significance threshold
[sigy, sigx] = find(pval<sig);

% colormap
load('helper_functions/redbluecmap');

N = size(rho,1);

imagesc(rho,[-1 +1]);
colormap(gca,flip(redbluecmap));
cbar = colorbar();
hold on
plot(sigx,sigy,'k*','MarkerSize',4);
hold off;
if ~exist('flag_labels','var')
    set(gca,'XTick',1:N,'YTick',1:N);
    xlabel('timeseries #');
    ylabel('timeseries #');
else
    set(gca,'XTick',[],'YTick',[]);
end
set(gca,'FontSize',14);

end