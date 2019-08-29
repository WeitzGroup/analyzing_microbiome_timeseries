
tmpdir = dir('data_mat/reg*');

% regression matrices - example at a single depth
fig = figure();
depthID = 1;
load(sprintf('data_mat/%s',tmpdir(depthID).name));
tmp_beta_label = {sprintf('Regression on X_{res}'),sprintf('Regression on X_{res}\nwith L1 norm')};
k = [2 4];
for kk = 1:length(k)
    subplot(1,2,kk);
    tmpim = imagesc(beta{k(kk)},[-1 +1]);
    tmpim.AlphaData = abs(beta{k(kk)})>betafrac_thresh;
    colormap(flip(redblue()));
    set(gca,'XTickLabel',[],'YTickLabel',[]);
    cbar = colorbar();
    cbar.Label.String = 'regression coefficient';
    title(tmp_beta_label{kk});
end
%fig.Position(3:4) = fig.Position(3:4).*[1 0.5];
fig.Position(3:4) = fig.Position(3:4).*[1.5 .75];
label_figure('A)',-0.2);
print_manuscript_fig('aloha3a');



% positive vs negative comparisons
fig = figure();
nrow = 2;
ncol = 4;
for depthID = 1:length(tmpdir)
	load(sprintf('data_mat/%s',tmpdir(depthID).name));
    subplot(nrow,ncol,depthID)
    bar(betafrac','stacked');
    set(gca,'XTickLabel',beta_label,'XTickLabelRotation',90);
    ylabel('fraction of pairs');
    title(sprintf('depth %dm',depth(depthID)));
    set(gca,'FontSize',12);
end

% set legend manually (......)
legend({'negative','positive'},'Location','EastOutside');
set(gca,'Position',[0.5570 0.2302 0.1419 0.2181]);

fig.Position(3:4) = fig.Position(3:4).*[1.75 1.25];
print('figures/regression/summary_posneg','-dpng');
label_figure('B)');
print_manuscript_fig('aloha3b');