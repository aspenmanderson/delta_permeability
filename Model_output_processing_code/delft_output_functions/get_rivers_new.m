function [CN, Tcl, L, Ncw, bw_centerline_noponds, channel_width_map] = get_rivers_new(input, shoreline_pixels, dx, len)

% Maps the river network and returns a logical matrix of where the rivers
% are located (uses RivaMap)
    %F. Isikdogan, A.C. Bovik, and P. Passalacqua. Automatic analysis of 
    %channel networks on remotely sensed images by singularity analysis,
    %IEEE Geosci. Remote Sens. Lett., in review.
    
    % Artificially fill the pixels around the shoreline
        % so there aren't artificial rivers along the shoreline
    for kk = 1:length(shoreline_pixels)
          for j = 1:length(shoreline_pixels{kk}(:,1))
                if shoreline_pixels{kk}(j,2) > 5
                    value = 5;
                elseif shoreline_pixels{kk}(j,2) > 4
                    value = 4;
                elseif shoreline_pixels{kk}(j,2) > 3
                    value = 3;
                elseif shoreline_pixels{kk}(j,2) > 2
                    value = 2;
                elseif shoreline_pixels{kk}(j,2) > 1
                    value = 1;
                end
                for i = 0:value
                    input(shoreline_pixels{kk}(j,2)-i,shoreline_pixels{kk}(j,1)) =  input(shoreline_pixels{kk}(j,2)-(value),shoreline_pixels{kk}(j,1));
                end
          end
    end
    
     % Artificially fill the pixels around the southern boundary
        % so there aren't artificial rivers along the shoreline
    [N, M] = size(input); %
    for kk = 1:M
        input(1,kk) = -5; %raise elevvation so rivers don't appear 
    end
   
    addpath("C:\Users\labuser10\OneDrive - sfu.ca\Working\SFU\Research\Modeling\delta_modeling_paper_codes\output_processing\delft_output_functions\RivaMap\"); %D3D

    % Run RivaMap
    [orientation, scaleMap, pos_psi, centerline] =  ModifiedSingularityIndex2D(input, 5, 1); %shown in test_code.m
    centerline  =  imfill(centerline);
    % this creates the filters that are needed to compute the multiscale
    % singularity index and applies the index to extract curvilinear structures
    % from the input. The singulatiy index function returns the overall
    % sinularoty index responce, width estimates, and channel orientation
    % for each pixel whether or not they are river centerlines

    % quiver plot (illustrates channel orientation)
    %figure;
    %imagesc(water_height); hold on
    [Nx, Ny] = size(input);
    xidx = 1:10:Nx;
    yidx = 1:10:Ny;
    [X,Y] = meshgrid(xidx,yidx);
    u = (pos_psi.*cos(orientation));
    u = u(xidx, yidx);
    v = (pos_psi.*sin(orientation));
    v = v(xidx, yidx);
    %quiver(Y',X',v,u, 'ShowArrowHead','on', 'LineWidth', 1, 'AutoScaleFactor', 2,'Color','k');

    % threshold centerlines
    centerline = mmscale(centerline);
    % extracts and theshold centerlines to delineate rivers
    level = graythresh(centerline);
    [gtT1r, gtT1c] = find(centerline > level*0.01); 
    gtT2 = centerline > level*0.1;
    bw_centerline = bwselect(gtT2, gtT1c, gtT1r, 8);

    %{
    figure;
    imagesc(bw_centerline); hold on
    colormap(gray);
    for k = 1:length(shoreline_pixels)
        plot(shoreline_pixels{k}(:,1), shoreline_pixels{k}(:,2), 'r','LineWidth', 2); hold on   
    end
    %}
  
    
    %Remove small connected components (e.g. ponds, tidal rivers?)
    %[Nx, Ny] = size(bw_centerline);
    CC = bwconncomp(bw_centerline);
    bw_centerline_noponds = false(Nx,Ny);
    minlength = len/10; %no reason given for min length in RivaMap

   
    for i = 1:length(CC.PixelIdxList)
        if(length(CC.PixelIdxList{i})>minlength)
            bw_centerline_noponds(CC.PixelIdxList{i}) = true;
        end
    end
    %{
    figure;
    imagesc(bw_centerline_noponds); hold on
    colormap(gray)
    %}
    
    % Regrow channels as raster
    [rowx,coly] = find(bw_centerline_noponds);
    Iregrown = bw_centerline_noponds;
    llen = scaleMap(bw_centerline_noponds) * 1.25;
    lthe = orientation(bw_centerline_noponds);

    x_off = -llen .* cos(lthe);
    y_off = llen .* sin(lthe);
    lines_all = [coly-x_off coly+x_off rowx-y_off rowx+y_off];

    for i = 1:length(lines_all)
        [xsn, ysn] = bresenham(lines_all(i,1), lines_all(i,2), lines_all(i,3), lines_all(i,4));
        linearInd = sub2ind(size(Iregrown), max(min(ysn,size(Iregrown,1)),1), max(min(xsn,size(Iregrown,2)),1));
        Iregrown(linearInd) = 1;
    end
    
    CC_noponds = bwconncomp(bw_centerline_noponds);

    % calculate width and length of the channels for each segment (Modified by Aspen)
    CN = 0; CW = []; CL = []; CL_pixels = []; channel_width = [];
    channel_width_map = zeros(Nx,Ny);
    clear rowx; clear coly;
    for i = 1:length(CC_noponds.PixelIdxList)

        % reset for next loop
        river = false(Nx,Ny);
        rowx = []; coly = [];

        % find each river segment 
        river(CC_noponds.PixelIdxList{i}) = true;
        [rowx,coly] = find(river);
        llen = scaleMap(river) * 1.25;
        lthe = orientation(river);

        % find the sides of the channel
        x_off = -llen .* cos(lthe);
        y_off = llen .* sin(lthe);
        lines = [coly-x_off coly+x_off rowx-y_off rowx+y_off];

        % calculate CN and channel width at outlet of CN (active channel does not need to intersect shoreline)
        [a,~] = size(rowx);
        for ii = 1:a
            channel_width(ii,i) = sqrt((lines(ii,2)-lines(ii,1))^2 + (lines(ii,4)-lines(ii,3))^2); %in m?
        end   
            
        if sum(channel_width(:,i)) > 0 && length(nonzeros(channel_width(:,i))) > 1 % making sur the river isn't just one cell
            CN = CN + 1;
            CW(CN) = median(nonzeros(channel_width(:,i)));
            CL_pixels(CN) = length(nonzeros(channel_width(:,i)));
            CL(CN) = max(max(bwdistgeodesic(river,coly(1),rowx(1))))*dx; % scale by cell size to put in [M]
            for ii = 1:a
                channel_width_map(rowx(ii),coly(ii)) = channel_width(ii,i);
            end 
        end        
    end  
  
    % incase there aren't any channels
    if CN == 0  
        CW(1) = 0;
        CL(1) = 0;
    end
    
    % sum total channel width divided by channel number to get an average
    % channel width for the delta
    Ncw = sum(CW)/CN;

    % find the maximum channel length (length of the main channel steam)
    L = max(CL);
    Tcl = sum(CL);

end

