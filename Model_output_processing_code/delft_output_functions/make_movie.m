function [outputArg1,outputArg2] = make_movie(inputArg1,inputArg2)
%UNTITLED3 Summary of this function goes here
%  Not finished
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Movie out output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%v = VideoWriter(strcat(ex_pth,model));
%v.FrameRate = 10; if length(lines_all) > 0 && area_condition > 20
    set(gcf,'position',[1,1,5*M_clip,5*N_clip])
    fig = imagesc(elev_clip); hold on
    colormap(gray);
    set(gca,'Xtick',1:50:401,'XTickLabel',{'1','10', '20', '30', '40','50','60','70','80'})
    set(gca,'Ytick',1:25:151,'YTickLabel',{'1','5','10','15','20','25','30'})
    xlabel('km'); ylabel('km');

    for k = 1:length(shoreline_pixels)
       h(1) = scatter(shoreline_pixels{k}(:,1), shoreline_pixels{k}(:,2),10,'filled','o','MarkerFaceColor',[1.00,1.00,0.00]); hold on 
       %C = sortrows(shoreline_pixels{k},1);
       %plot(C(:,1), C(:,2),'Color',[0.25,0.25,0.25],'Linewidth',1.5); hold on
    end

    h(2) = plot(x_clino, clino_profile,'Linewidth',1.5,'Color',[0.00,1.00,0.40]); hold on %[0.40,0.75,0.46]

    for i = 1:length(lines_all)
       h(3) = line(lines_all(i,1:2),lines_all(i,3:4),'Color','b','LineWidth', 0.75); hold on % [0.00,0.40,0.80]
    end

    set(gca,'YDir','normal')
    g = colorbar; %caxis([0 1]);
    set(get(g,'label'),'string','Elevation (masl)','Rotation',90.0);

    title(['Time elapsed in years = ' num2str(round(morft{step}*morfac/360))]) 
    legend(h([1, 2, 3]),{'\color{white} shoreline','\color{white} subaqueous clinoform','\color{white} river nerwork'})
    set(legend,'color','none'); legend boxoff  

    % annotation
    dim = [0.15 0.6 0.3 0.3];
    str = {['AD: ' num2str(round(AD(ii),1))] , ['CN: ' num2str(round(CN(ii),0))],['Tcw: ' num2str(round(Tcw(ii),1))],['Dsh: ' num2str(round(Dsh(ii),1))]};
    t = annotation('textbox',dim,'String',str,'FitBoxToText','on');
    t.FontSize = 16;
    t.Color = 'w';
    t.FontName = 'Times';
    t.EdgeColor = 'none';

    frame = getframe(gcf);
    %frame = imresize(frame,[1924 794]);
    writeVideo(v,frame);

    my_cmap2 = load(strcat(cp_pth,'aa_delft1.mat'));

    set(gcf,'position',[1,1,5*M_clip,5*N_clip])
    fig = imagesc(elev_clip); hold on
    colormap(gca,(my_cmap2.cmap));
    set(gca,'Xtick',1:50:401,'XTickLabel',{'1','10', '20', '30', '40','50','60','70','80'})
    set(gca,'Ytick',1:25:151,'YTickLabel',{'1','5','10','15','20','25'})
    xlabel('km'); ylabel('km');

    set(gca,'YDir','normal')
    h = colorbar; %caxis([0 1]);
    set(get(h,'label'),'string','Elevation (masl)','Rotation',90.0);

    title(['Time elapsed in years = ' num2str(round(morft{step}*morfac/360))]) 

    frame = getframe(gcf);
    writeVideo(v,frame);

    time_step(ii) = step;
    morpho_time(ii) = morft{step}*morfac/360;
    time_condition = morpho_time(ii);

    ii = ii + 1;
    %delete(findall(gcf,'type','annotation'))
%}
end

