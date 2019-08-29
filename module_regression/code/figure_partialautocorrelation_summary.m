
clear;
tmpdir = dir('data_mat/pac*');

nrow = 2;
ncol = 4;

% average pac for each depth 
fig = figure();
for depthID = 1:length(tmpdir)
    load(sprintf('data_mat/%s',tmpdir(depthID).name));  
    subplot(nrow,ncol,depthID);
    tmppac = apply_threshold(pac,pac_thresh);
    boxplot(tmppac',lag);
    ylim([-1 1]);
    xlabel('lag');
    ylabel('PAC');
    title(sprintf('depth %dm',depth(depthID)));
    set(gca,'FontSize',12);
end
fig.Position(3) = fig.Position(3)*1.75;
print('figures/partialautocorrelation/summary_value','-dpng');
label_figure('B)');
print_manuscript_fig('aloha2b');

% how many contigs show autocorrelation? (by depth)
fig = figure();
for depthID = 1:length(tmpdir)
    load(sprintf('data_mat/%s',tmpdir(depthID).name));  
    subplot(nrow,ncol,depthID);
    tmpfreq = abs(pac)<=1 & abs(pac)>pac_thresh;
    tmpfreq = sum(tmpfreq,2)/length(contigID);
    bar(lag,tmpfreq);
    ylim([0 1]);
    xlabel('lag');
    ylabel('fraction of contigs');
    title(sprintf('depth %dm',depth(depthID)));
    set(gca,'FontSize',12);
end
fig.Position(3) = fig.Position(3)*1.75;
print('figures/partialautocorrelation/summary_fraction','-dpng');
label_figure('A)');
print_manuscript_fig('aloha2a');

% omit pacs>1 and pacs<thresh
function pac = apply_threshold(pac,thresh)
    pac(abs(pac)>1) = nan;
    pac(abs(pac)<thresh) = nan;
end