
tmpdir = dir('data_mat/reg*');

for depthID = 1:length(tmpdir)
    
    % load data
    load(sprintf('data_mat/%s',tmpdir(depthID).name));
    fig = figure();
    nrow = 1;
    ncol = length(beta);
    
    % plot each regression analysis result (contained in the cell beta)
    for k = 1:length(beta)
        subplot(nrow,ncol,k);
        tmpim = imagesc(beta{k},[-1 +1]);
        tmpim.AlphaData = abs(beta{k})>betafrac_thresh; % zero interactions
        colormap(flip(redblue()));
        colorbar();
        set(gca,'XTickLabel',[],'YTickLabel',[]);
        title(beta_label{k});
    end
    
    % save figure & close
    fig.Position(3:4) = fig.Position(3:4).*[2 .5];
    print(sprintf('figures/regression/depth%d',depthID),'-dpng');
    close;
    
end

