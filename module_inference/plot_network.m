function plot_network(X,zero_thresh,fs)

im = imagesc(X);
im.AlphaData = X>zero_thresh;
colormap(gca,'jet');
colorbar();
set(gca,'ytick',[],'xtick',[]);
set(gca,'FontSize',fs);
%title('$\tilde{M}_{rec}$','Interpreter','latex','FontSize',fs+2);
axis square;
xlabel('viruses');
ylabel('hosts');

end