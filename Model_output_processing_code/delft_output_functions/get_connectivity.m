function [C, NumC, SizeC, percent_hp_and_cn] = get_connectivity(input, threshold, channel_network_map, AD, dx, elev_clip)

% Maps the high permeabable areas and calculates a connectivity value
    % Steel et al. 2022

    % convert perm map to binary image based on threshold
    [M,N] = size(input);
    bwimag = zeros(M,N);
    permimag = nan(M,N);
    
    for i = 1:M
        for j = 1:N
            if input(i,j) >= threshold
                bwimag(i,j) = 1;
                permimag(i,j) = input(i,j);
            end
        end
    end
    
    %{
    figure
    imagesc(bwimag); hold on
    colormap(gray)
    for k = 1:length(shoreline_pixels)
        plot(shoreline_pixels{k}(:,1), shoreline_pixels{k}(:,2), 'r','LineWidth', 2); hold on   
    end
    %}
    %% 
    %{
    hf = figure;
    h1 = axes;
    p1 = imagesc(elev_clip);%,'AlphaData',imAlpha); hold on
    colormap(h1,gray(225));
    set(h1,'YDir','normal')
    alpha(h1,0.5)
    hold on
   % colorbar
    %set(h1,'Xtick',[]); set(h1,'Ytick',[])
    set(gca,'Xtick',[])
    set(gca,'Ytick',[])
    %
   
   
    h2 = axes;  
    imAlpha = ones(size(permimag));
    imAlpha(isnan(permimag)) = 0;
    p2 = imagesc(permimag,'AlphaData',imAlpha); 
    set(h2,'color', 'none','visible','off')%'color',[1 1 1]); % sets background as white
    
    cp_pth = ("C:\Program Files\Deltares\Delft3D 4.03.01\win64\quickplot\bin\colormaps\"); %colorplot location
    my_cmap = load(strcat(cp_pth,'perm3_cmap.mat'));
    colormap(h2,my_cmap.cm);
    caxis([5E-1 5E2])
    caxis([7E-2 1E1])
    set(h2,'ColorScale','log')
    set(h2,'YDir','normal')
    set(gca,'Xtick',[])
    set(gca,'Ytick',[])
    linkaxes([h1 h2])
   
    set(h2,'Xtick',[]); set(h2,'Ytick',[])
    set(gcf,'position',[100 100 800 450])
    %}
%%
   
    % caxis([3E-9 8E-8]) %for hydrualic conductivity
    set(gca,'YDir','normal')
    %h = colorbar; %caxis([0 1]);
   % set(get(h,'label'),'string','\kappa(Darcy)','Rotation',90.0);
       
    % Find connected components
    CC = bwconncomp(bwimag);
    
    % Calculate Connectivity
    if CC.NumObjects > 0 % make sure there are high perm bodies
         % number of pixels in the largest body 
         SAb = 0;
         Ab = length(CC.PixelIdxList{1});
         for i = 1:length(CC.PixelIdxList)
            if length(CC.PixelIdxList{i}) > Ab
                Ab = length(CC.PixelIdxList{i});
            end
         % sum of high perm pixels 
            SAb = SAb + length(CC.PixelIdxList{i});
         end
         C = Ab/SAb;
         SizeC_m2 = Ab*dx*dx; % number of land pixels times the size of the pixel [M^2]
         SizeC = ((SizeC_m2*1e-6)/AD)*100; 
   else
        C = NaN;
        SizeC = 0;
   end
   
   NumC = CC.NumObjects;
   
   % Define percentage of the high permeability areas overlap with the
   % channel network
   compositimag = zeros(M,N);
   for i = 1:M
       for j = 1:N
           if bwimag(i,j) == 1 && channel_network_map(i,j) == 1
                compositimag(i,j) = 1;
            end
        end
   end
   
    %{
    figure
    imagesc(compositimag); hold on
    colormap(gray)
    for k = 1:length(shoreline_pixels)
        plot(shoreline_pixels{k}(:,1), shoreline_pixels{k}(:,2), 'r','LineWidth', 2); hold on   
    end
    %}
   
   Scn = sum(sum(channel_network_map == 1));
   Shp_cn = sum(sum(compositimag == 1));
   
   percent_hp_and_cn = NaN;
   if Scn > 0 
       percent_hp_and_cn = (Shp_cn/Scn)*100;
   end
   
end

