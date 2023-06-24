function [range,sill,gamma] = get_semivarience(input_mask,dx,dy,model)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calculate varience, covarience, and semivarigram of permeability 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ex_pth = ("I:\Aspen_Delft3D_local\Summer_2020_trials\output_CC_try6_REDO\figures_semivarience\");
    cp_pth = ("C:\Program Files\Deltares\Delft3D 4.03.01\win64\quickplot\bin\colormaps\"); %colorplot location
    addpath ('C:\Users\labuser10\OneDrive - sfu.ca\Working\SFU\Research\Modeling\delta_modeling_paper_codes\output_processing\delft_output_functions') %variogram function
    addpath ('C:\Users\labuser10\OneDrive - sfu.ca\Working\SFU\Research\Modeling\delta_modeling_paper_codes\output_processing\delft_output_functions\variogramfit') %variogram characteristics function
    
    % flip and transpose matrix so azimuth angles are from the horizontal,
    % not from 90 degrees 
    %input_flip = flipud(input_mask);
   % input_T = input_flip.';
    input_T = flipud(input_mask);
    
    % calculate position vector(location of data and clip the edges of the variable 
    [X,Y] = meshgrid(1:length(input_mask(1,:)),1:length(input_mask(:,1)));
    X = X*dx;
    Y = Y*dy;
    %perm_mask_T = perm_mask.';
    val = []; pos = [];
    k = 1;
    for i = 1:length(X(:,1))
        for j = 1:length(X(1,:))
            pos(k,1) = X(i,j);  % pos = [ndata,ndims]
            pos(k,2) = Y(i,j); 
            val(k,1) = input_T(i,j);
            k = k + 1;
        end
    end
    nbins = floor(length(Y)/10); %20 works and runs a lot faster
    
    
    %figure
    %my_cmap = load(strcat(cp_pth,'directional_vario.mat'));
    
    % experimental semivariogram
    gamma = variogram(pos,val,'nrbins',nbins,'anisotropy',true,'thetastep',30,'plotit',false,'subsample',10000);
    % fit theorestical variogram to experimental
    
     theta = rad2deg(gamma.theta);
     nugget = []; range = []; sill = [];
     for i  = 1:length(gamma.val(1,1:end))-1
        [a,c,n,S] = variogramfit(gamma.distance, gamma.val(:,i));
        %nugget(i) = S.nugget;
        range(i) = S.range;
        sill(i) = S.sill;
     end
    
    %{
    fig = gcf;
    set(fig,'Position',[100 100 800 400])
    view(0,90)
    xlabel('Lag Distance (m): y-direction')
    ylabel('Lag Distance (m): x-direction')
    h = colorbar; %caxis([0 1]);
    ylim([0 20000])
    xlim([-25000 25000])
    %colormap(gca,my_cmap.cm);
    %colormap(gca,turbo);
    set(get(h,'label'),'string','\gamma (h)','Rotation',90.0);
    %if mod_num ~= 3
        %caxis([1E-14 2.3E-13])
    %else
        caxis([4.5E-15 2.3E-13])
    %end
    set(gca,'colorscale','log')

   % saveas(fig,strcat(ex_pth,model,'_dirvario.eps'));
   % saveas(fig,strcat(ex_pth,model,'_dirvario.jpg'));
   % saveas(fig,strcat(ex_pth,model,'_dirvario.pdf'));
   % saveas(fig,strcat(ex_pth,model,'_dirvario.fig'));
    saveas(fig,strcat(ex_pth,model,'_dirvario.png'));
    %}
    line_colors(1:3,:) = [0.18, 0.85, 0.52; 0.31, 0.31, 0.31; 0.93, 0.88, 0.12];
    line_colors(4:6,:) = [0.18, 0.22, 0.69; 0.35, 0.56, 0.35; 0.06, 0.65, 0.88];
    line_colors(7,:) = [0.18, 0.85, 0.52];

    figure
    Markers = {'x'; 'square'; 'o'; 'diamond'; 'pentagram'; 'hexagram'; 'x'};
    for i  = 1:length(gamma.val(1,1:end))
        plot(gamma.distance,gamma.val(:,i),'MarkerSize',5,'Marker',Markers{i},'LineWidth',1.5,'Color',line_colors(i,:)); hold on
    end
    fig = gcf
    set(fig,'Position',[100 100 800 350])
   
    xlim([min(gamma.distance) max(gamma.distance)])
    %if mod_num == 3
    %    legend('0','30','60','90','120','150','180','location','northwest')
    %end
    %legend boxoff
    xlabel('Lag Distance (m)')
    ylabel('\gamma (h)')
    %if mod_num ~= 3
    %    ylim([0 2.3E-13])
    %else
    %    ylim([0 2.3E-14])
    %end
    %ylim([7E-16 2.3E-13])
    set(gca, 'YScale', 'log')
    
   % saveas(fig,strcat(ex_pth,model,'_semivario.eps'));
   % saveas(fig,strcat(ex_pth,model,'_semivario.jpg'));
    saveas(fig,strcat(ex_pth,model,'_semivario.pdf'));
   % saveas(fig,strcat(ex_pth,model,'_semivario.fig'));
   % saveas(fig,strcat(ex_pth,model,'_semivario.png'));
end

