function label_ax(str,xpos)
% add figure label

if ~exist('xpos','var') || isempty(xpos)
    xpos = -0.15;
end

text(gca,xpos,1.05,str,'Units','normalized','FontSize',8,'FontWeight','bold');

end