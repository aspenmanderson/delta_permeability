function [CN, Tcw, L, lines_all, bw_centerline_noponds] = get_rivers(input, shoreline_pixels, dx)

% Maps the river network and returns a logical matrix of where the rivers
% are located (uses RivaMap)
    %F. Isikdogan, A.C. Bovik, and P. Passalacqua. Automatic analysis of 
    %channel networks on remotely sensed images by singularity analysis,
    %IEEE Geosci. Remote Sens. Lett., in review.
    
    
    % Artificially fill the pixels around the shoreline
        % so there aren't artificial rivers along the shoreline
    [N, M] = size(input); 
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
    %}
   
    addpath("C:\Users\labuser10\OneDrive - sfu.ca\Working\SFU\Research\Modeling\delta_modeling_paper_codes\output_processing\delft_output_functions\RivaMap\"); %D3D

    % Run RivaMap
    
    [orientation, scaleMap, pos_psi, centerline] =  ModifiedSingularityIndex2D(input, 5, 1); %shown in test_code.m
    centerline =  imfill(centerline);
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
    minlength = 10; %no reason given for min length in RivaMap

    %figure;
    for i = 1:length(CC.PixelIdxList)
        if(length(CC.PixelIdxList{i})>minlength)
            bw_centerline_noponds(CC.PixelIdxList{i}) = true;
        end
    end
    %imagesc(bw_centerline_noponds); hold on
    %colormap(gray)
%
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
    %
    CC_noponds = bwconncomp(bw_centerline_noponds);

    %{
    figure;
    imagesc(water_height); hold on
    %clear line;
    for i = 1:length(lines)
        line(lines(i,1:2), lines(i,3:4));
    end
    colormap(gray)
    %}

    % calculate width and length of the channels for each segment (Modified by Aspen)
    channel_width = []; CN = 0; CW = []; chan_length = 0;
    channel_width_map = zeros(Nx,Ny);
    clear rowx; clear coly;
    %active_river = 1;
    %CN = length(CC_noponds.PixelIdxList);
    for i = 1:length(CC_noponds.PixelIdxList)

        % reset for next loop
        river = false(Nx,Ny);
        channel_width = []; rowx = []; coly = [];

        % find each river segment 
        river(CC_noponds.PixelIdxList{i}) = true;
        [rowx,coly] = find(river);
        llen = scaleMap(river) * 1.25;
        lthe = orientation(river);

        % find the sides of the channel
        x_off = -llen .* cos(lthe);
        y_off = llen .* sin(lthe);
        lines = [coly-x_off coly+x_off rowx-y_off rowx+y_off];

        % calculate channel width
        [a,~] = size(rowx);
        for ii = 1:a
            channel_width(ii) = sqrt((lines(ii,2)-lines(ii,1))^2 + (lines(ii,4)-lines(ii,3))^2); %in m?
            channel_width_map(rowx(ii),coly(ii)) = channel_width(1,ii);
        end
        
        % old way of finding CN only at shoreline
        j = 1; %starts at 2 to avoid checking the cells along the southern and west boundaries
        for kk = 1:length(shoreline_pixels)
            while j <= length(shoreline_pixels{kk}(:,1))
                
                if river(shoreline_pixels{kk}(j,2),shoreline_pixels{kk}(j,1)) == 1 %shoreline pixel
                     cw = channel_width_map(shoreline_pixels{kk}(j,2),shoreline_pixels{kk}(j,1));
                     if cw > 0
                         CN = CN + 1; % add 1 for each active channel   
                         CW(CN) = cw;
                         j = j + 2*ceil(CW(CN)); %jump ahead one channel width 
                     else
                         j = j + 1;
                     end
                     continue
                     
                %south of shoreline pixel    
                elseif river(shoreline_pixels{kk}(j,2)-1,shoreline_pixels{kk}(j,1)) == 1 && shoreline_pixels{kk}(j,1)-1 > 0 
                     cw = channel_width_map(shoreline_pixels{kk}(j,2)-1,shoreline_pixels{kk}(j,1));
                     if cw > 0
                         CN = CN + 1; % add 1 for each active channel   
                         CW(CN) = cw;
                         j = j + 2*ceil(CW(CN)); %jump ahead one channel width 
                     else
                         j = j + 1;
                     end
                     continue
                     
                %west of shoreline pixel     
                elseif river(shoreline_pixels{kk}(j,2),shoreline_pixels{kk}(j,1)-1) == 1 && shoreline_pixels{kk}(j,1)-1 > 0
                     cw = channel_width_map(shoreline_pixels{kk}(j,2),shoreline_pixels{kk}(j,1)-1);
                     if cw > 0
                         CN = CN + 1; % add 1 for each active channel   
                         CW(CN) = cw;
                         j = j + 2*ceil(CW(CN)); %jump ahead one channel width 
                     else
                         j = j + 1;
                     end
                     continue
                
                %southwest of shoreline pixel
                elseif river(shoreline_pixels{kk}(j,2)-1,shoreline_pixels{kk}(j,1)-1) == 1 && shoreline_pixels{kk}(j,1)-1 > 0
                     cw = channel_width_map(shoreline_pixels{kk}(j,2)-1,shoreline_pixels{kk}(j,1)-1);
                     if cw > 0
                         CN = CN + 1; % add 1 for each active channel   
                         CW(CN) = cw;
                         j = j + 2*ceil(CW(CN)); %jump ahead one channel width 
                     else
                         j = j + 1;
                     end
                     continue
                     
                    if shoreline_pixels{kk}(j,1)+1 <= length(input(1,:)) % check to make sure the river pixel is not out of bounds(larger than the model) 
                        elseif river(shoreline_pixels{kk}(j,2),shoreline_pixels{kk}(j,1)+1) == 1; %east of shoreline pixel
                             cw = channel_width_map(shoreline_pixels{kk}(j,2),shoreline_pixels{kk}(j,1)+1);
                             if cw > 0
                                 CN = CN + 1; % add 1 for each active channel   
                                 CW(CN) = cw;
                                 j = j + 2*ceil(CW(CN)); %jump ahead one channel width 
                             else
                                 j = j + 1;
                             end
                             continue 
                        elseif river(shoreline_pixels{kk}(j,2)-1,shoreline_pixels{kk}(j,1)+1) == 1 %southeast of shoreline pixel
                             cw = channel_width_map(shoreline_pixels{kk}(j,2)-1,shoreline_pixels{kk}(j,1)+1);
                             if cw > 0
                                 CN = CN + 1; % add 1 for each active channel   
                                 CW(CN) = cw;
                                 j = j + 2*ceil(CW(CN)); %jump ahead one channel width 
                             else
                                 j = j + 1;
                             end
                             continue
                    end
                    continue

                    if shoreline_pixels{kk}(j,2)+1 <= length(input(:,1)) % check to make sure the river pixel is not out of bounds(larger than the model)
                        elseif river(shoreline_pixels{kk}(j,2)+1,shoreline_pixels{kk}(j,1)) == 1 %north of shoreline pixel
                             cw = channel_width_map(shoreline_pixels{kk}(j,2)+1,shoreline_pixels{kk}(j,1));
                             if cw > 0
                                 CN = CN + 1; % add 1 for each active channel   
                                 CW(CN) = cw;
                                 j = j + 2*ceil(CW(CN)); %jump ahead one channel width 
                             else
                                 j = j + 1;
                             end
                             continue
                        elseif river(shoreline_pixels{kk}(j,2)+1,shoreline_pixels{kk}(j,1)+1) == 1 %northeast of shoreline pixel
                             cw = channel_width_map(shoreline_pixels{kk}(j,2)+1,shoreline_pixels{kk}(j,1)+1);
                             if cw > 0
                                 CN = CN + 1; % add 1 for each active channel   
                                 CW(CN) = cw;
                                 j = j + 2*ceil(CW(CN)); %jump ahead one channel width 
                             else
                                 j = j + 1;
                             end
                             continue    
                        elseif river(shoreline_pixels{kk}(j,2)+1,shoreline_pixels{kk}(j,1)-1) == 1 %northwest of shoreline pixel
                             cw = channel_width_map(shoreline_pixels{kk}(j,2)+1,shoreline_pixels{kk}(j,1)-1);
                             if cw > 0
                                 CN = CN + 1; % add 1 for each active channel   
                                 CW(CN) = cw;
                                 j = j + 2*ceil(CW(CN)); %jump ahead one channel width 
                             else
                                 j = j + 1;
                             end
                             continue
                    end
                end
                j = j + 1;
            end
        end
        
        % incase there aren't any channels
        if CN == 0  
            CW(1) = 0;
        end
                        
        % calculate number of cells in each channel
        chan_length(i) = max(max(bwdistgeodesic(river,coly(1),rowx(1))))*dx; % scale by cell size to put in [M]

    end

    % sum total channel width
    Tcw = sum(CW);

    % find the maximum channel length (length of the main channel steam)
    L = max(chan_length);

end

