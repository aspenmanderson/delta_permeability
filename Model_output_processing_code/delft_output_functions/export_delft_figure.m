function [] = export_delft_figure(M_clip,N_clip,elev_clip,morft,ex_pth,model,label)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    cp_pth = ("C:\Program Files\Deltares\Delft3D 4.03.01\win64\quickplot\bin\colormaps\"); %colorplot location
    my_cma6 = load(strcat(cp_pth,'aa_delft7.mat'));
    my_cmap = load(strcat(cp_pth,'aa_delft3.mat'));
    
    figure
    
    set(gcf,'position',[100 100 800 450]) 
    %set(gcf,'position',[1,1,3*M_clip,5*N_clip])
    fig = imagesc(elev_clip); hold on
    colormap(gca,(my_cma6.cmap3));
    %colormap(gca,my_cmap.cm);
    if label == 1
        set(gca,'Xtick',1:50:M_clip,'XTickLabel',{'1','10', '20', '30', '40','50','60','70','80'})
        set(gca,'Ytick',1:50:N_clip,'YTickLabel',{'1','10', '20', '30', '40','50','60','70','80'})
        xlabel('km'); ylabel('km');
    elseif label == 0 
        set(gca,'Xtick',[])
        set(gca,'Ytick',[])
    end
    

    set(gca,'YDir','normal')
    h = colorbar;
    %caxis([-35 15]);
    caxis([-15 5]);
    %caxis([-25 12])
    set(get(h,'label'),'string','Elevation (masl)','Rotation',90.0);

%    title(['Time elapsed in years = ' num2str(round(morft{step}/360))]) 
    %saveas(fig,strcat(ex_pth,'figures\',model,'.eps'));
    saveas(fig,strcat(ex_pth,'figures\',model,'.jpg'));
    %saveas(fig,strcat(ex_pth,'figures\',model,'.png'));
    saveas(fig,strcat(ex_pth,'figures\',model,'.pdf'));
    %saveas(fig,strcat(ex_pth,'figures\',model,'.fig'));

end

