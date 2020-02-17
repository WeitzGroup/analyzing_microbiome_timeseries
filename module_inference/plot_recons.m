function plot_recons(t, H, V, h, W, Mtilde, Mtilde_hat)

    clf
    fs = 10;
    lw = 1.5;
    fig = gcf;
    fig.Units = 'inches';
    fig.Position(3:4) = [6.5 7];
    fig.Position(1:2) = [7 2];
    
    label_ax = @(str) text(-0.4,1.6,str,'FontWeight','bold','FontSize',fs+2,'Units','inches');
    
    ax(1) = subplot(3,2,1);
    semilogy(t,H,'LineWidth',lw);
    xlim([t(1) t(end)]);
    tmpy = log10(ylim);
    tmpy(1) = floor(tmpy(1));
    tmpy(2) = ceil(tmpy(2));
    ylim(10.^tmpy);
    xlabel('time (hrs)');
    ylabel('density (1/mL)');
    title('host dynamics','FontWeight','normal');
    set(gca,'FontSize',fs);
    label_ax('a)');
    
    ax(2) = subplot(3,2,2);
    semilogy(t,V,'LineWidth',lw);
    xlim([t(1) t(end)]); 
    tmpy = log10(ylim);
    tmpy(1) = floor(tmpy(1));
    tmpy(2) = ceil(tmpy(2));
    ylim(10.^tmpy);
    xlabel('time (hrs)');
    %ylabel('density (1/mL)');
    title('virus dynamics','FontWeight','normal');
    set(gca,'FontSize',fs);
        
    ax(3) = subplot(3,2,3);
    imagesc(h);
    colorbar();
    set(gca,'ytick',[],'xtick',[])
    ylabel('hosts');
    set(gca,'FontSize',fs);
    title('H','Interpreter','latex','FontSize',fs+2);
    label_ax('b)');
    
    ax(4) = subplot(3,2,4);
    imagesc(W);
    colorbar();
    ylabel('viruses');
    colormap(gca,'hot');
    set(gca,'ytick',[],'xtick',[])
    set(gca,'FontSize',fs);
    title('W','Interpreter','latex','FontSize',fs+2);
    
    ax(5) = subplot(3,2,5);
    im = imagesc(Mtilde);
    im.AlphaData = Mtilde~=0;
    colormap(gca,'jet');
    colorbar();
    set(gca,'ytick',[],'xtick',[])
    set(gca,'FontSize',fs);
    title('$\tilde{M}$','Interpreter','latex','FontSize',fs+2);
    %axis square;
    label_ax('c)');
    
    zero_thresh = 1e-7;    % set threshold for 'zero'
    ax(6) = subplot(3,2,6);
    im = imagesc(Mtilde_hat);
    im.AlphaData = Mtilde_hat>zero_thresh;
    colormap(gca,'jet');
    colorbar();
    set(gca,'ytick',[],'xtick',[]);
    set(gca,'FontSize',fs);
    title('$\tilde{M}_{rec}$','Interpreter','latex','FontSize',fs+2);
    %axis square;
    tmpstr = sprintf('relative error = %.2g',inference_error(Mtilde_hat,Mtilde));
    text(1,-0.1,tmpstr,'Units','normalized','HorizontalAlignment','right','FontSize',fs);
    
    % adjust subplot positions
    ax(3).Position(2) = ax(3).Position(2)-0.04;
    ax(4).Position(2) = ax(4).Position(2)-0.04;
    ax(5).Position(2) = ax(5).Position(2)-0.04;
    ax(6).Position(2) = ax(6).Position(2)-0.04;   
    %ax(5).Position(1) = ax(5).Position(1)-0.04;
    %ax(6).Position(1) = ax(6).Position(1)-0.04;
end
        
