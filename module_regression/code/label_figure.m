function label_figure(str,xpos)
% add figure label

if ~exist('xpos','var') || isempty(xpos)
    xpos = -0.4;
end

tmpfig = gcf;
tmpax = tmpfig.Children(end);
text(tmpax,xpos,1.05,str,'Units','normalized','FontSize',16,'FontWeight','bold');

end