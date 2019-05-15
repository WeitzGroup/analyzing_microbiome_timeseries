
% all timeseries in one figure

mode_flag = 0; % area plots
%mode_flag = 1; % bar plots

tmpdir = dir('data_mat/ts*');
fig = figure();
nrow = 2;
ncol = 4;


for k = 1:length(tmpdir)
    load(sprintf('data_mat/%s',tmpdir(k).name));

    subplot(nrow,ncol,k);
    
	% get rid of nans so that area plots do not have gaps
    if ~exist('nanID','var')
        nanID = sum(isnan(X0))==length(contigID);
    end
    tmpt = t(~nanID);
    tmpX = X0(:,~nanID)';

    % sort so that most abundant are on bottom
    [~,tmpID] = sort(sum(tmpX,'omitnan'),'descend');
    tmpX = tmpX(:,tmpID);

    % plot abundance + sampling info
    if mode_flag==0
        area(tmpt,tmpX);
        hold on; % mark missing samples
        plot([t(nanID) t(nanID)]',repmat([0 1],[length(t(nanID)) 1])','w--','LineWidth',1.5);
        hold off;
    elseif mode_flag==1
        bar(tmpt,tmpX,1,'stacked');
    end

    % labels
    title(sprintf('depth %dm',depth(k)));
    ylabel('relative abundance');
    ylim([0 1]);
    xlim([tmpt(1) tmpt(end)]);
    set(gca,'XTick',t,'XTickLabel',string(t,'MMM-d-yy'),'XTickLabelRotation',90);
    set(gca,'FontSize',16);

end
%clear k tmp*;

% print
fig.Position(3:4) = fig.Position(3:4).*[3 1.5];
print('figures/timeseries/summary','-dpng');
print_manuscript_fig('aloha1');
