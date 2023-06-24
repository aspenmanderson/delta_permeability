function [] = export_perm_figure(M_clip,N_clip,perm_mask,morft,ex_pth,model,label)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    cp_pth = ("C:\Program Files\Deltares\Delft3D 4.03.01\win64\quickplot\bin\colormaps\"); %colorplot location
    my_cmap = load(strcat(cp_pth,'perm3_cmap.mat'));
    
    figure
    set(gcf,'position',[100 100 800 450]) 
    %set(gcf,'position',[1,1,5*M_clip,5*N_clip])
    imAlpha = ones(size(perm_mask));
    imAlpha(isnan(perm_mask)) = 0;
    fig = imagesc(perm_mask,'AlphaData',imAlpha); hold on
    set(gca,'color',[1 1 1]); % sets background as white
    colormap(gca,my_cmap.cm);
    set(gca,'ColorScale','log')
    set(gca,'Xtick',[])
    set(gca,'Ytick',[])
    %xlabel('km'); ylabel('km');
   % xlim([50 350]); ylim([0 125])
    
    
    if label == 1
        set(gca,'Xtick',1:50:M_clip,'XTickLabel',{'1','10', '20', '30', '40','50','60','70','80'})
        set(gca,'Ytick',1:50:N_clip,'YTickLabel',{'1','10', '20', '30', '40','50','60','70','80'})
        xlabel('km'); ylabel('km');
    elseif label == 0 
        set(gca,'Xtick',[])
        set(gca,'Ytick',[])
    end
    
   % caxis([5E-1 5E2])
    set(gca,'YDir','normal')
    h = colorbar; %caxis([0 1]);
    set(get(h,'label'),'string','\kappa(Darcy)','Rotation',90.0);

    title(['Time elapsed in years = ' num2str(round(morft/360))]) 
    %saveas(fig,strcat(ex_pth,'figures\',model,'.eps'));
    saveas(fig,strcat(ex_pth,'figures_perm\',model,'.jpg'));
    saveas(fig,strcat(ex_pth,'figures_perm\',model,'.png'));
    saveas(fig,strcat(ex_pth,'figures_perm\',model,'.pdf'));
    %saveas(fig,strcat(ex_pth,'figures\',model,'.fig'));

end

