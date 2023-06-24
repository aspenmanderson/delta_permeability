function [mask, shoreline_pixels] = get_shoreline(Water, threshold)

% Delineating the shoreline based on the Open Angle Method (Shaw et. al 2008)
    %Theory: a point that sees 180 degrees of open ocean is part of the
    %shoreline but a point that sees less than 180 degrees of open ocean is not
   
    warning ('off','all');
    
    % STEP 1: create a binary map with a set of water points (W) and land points (L)
        % (P) are all grid cells in the domain 
    [N, M] = size(Water); % resize grid without shore
    [X,Y] = meshgrid(1:1:M,1:1:N);
    
    %{
    Water = zeros([N,M]);
    for i = 1:N
        for j = 1:M
            if elev_clip(i,j) - wl_clip(i,j) < 0
               Water(i,j) = 1;
            end
        end
    end
    %}
    
    
%{
        

    P = [];
    Land = zeros([N,M]);
    L = [];
    Water = zeros([N,M]);
    W = [];
    for i = 1:N
        for j = 1:M
            P = [P; X(i,j) Y(i,j)];
            if elev_clip(i,j) - wl_clip(i,j) > 0
               Land(i,j) = 1;
               L = [L; X(i,j) Y(i,j)];
            else
               Water(i,j) = 1;
               W = [W; X(i,j), Y(i,j)];
            end
        end
    end
%}
   
    %{
    % Define the set of points that represent the land-water interface (I)
    Interface = NaN([N,M]);
    I = [];
    for i = 3:N-2 %indexing 
        for j = 3:M-2
            if Land(i,j) == 1
                if Land(i+1,j) == 0 || Land(i,j+1) == 0 || Land(i-1,j) == 0 || Land(i,j-1) == 0
                   Interface(i,j) = 1;
                   I = [I; X(i,j) Y(i,j)];
                elseif Land(i+1,j-1) == 0 || Land(i-1,j+1) == 0 || Land(i+1,j+1) == 0 || Land(i-1,j-1) == 0
                   Interface(i,j) = 1;
                   I = [I; X(i,j) Y(i,j)];
                end
            end
        end
    end

    %{
    % check what land-water interface looks like
    surface = Interface;
    my_cmap = load(strcat(cp_pth,'aa_slainity.mat'));
    figure; pcolor(surface); hold on
    shading(gca,'interp')
    colormap(gca,flipud(my_cmap.cm));
    plot(I(:,1),I(:,2),'.')
    %}

    
    % STEP 2: Find the convex hull of the land points (Lc)
    [idx, area] = convhull(L); %indicies of where the hull exists in L

    %{
    % check indicies of convex hull
    plot(L(idx,1),L(idx,2))
    plot(L(idx,1),L(idx,2),'*')
    %}

    % create a polygon for a convex hull (Lc)
    Lc_poly = polyshape(L(idx,:));
    idx = [];
    idx = isinterior(Lc_poly,P);
    Lc = [P(idx,1) P(idx,2)];
    
    %{
    % check what convex hull looks like
    surface = Land;
    my_cmap = load(strcat(cp_pth,'aa_slainity.mat'));
    figure; pcolor(surface); hold on
    shading(gca,'interp')
    colormap(gca,flipud(my_cmap.cm));
    plot(Lc_poly)
    %}

    % STEP 3: Define points of interst (W-Wo (Or faster computation, Lc - Land) and I)

    % Find Lc - L
    k = 1; j = 1; POI = [];
    for i = 1:length(Lc)
        if [Lc(i,1),Lc(i,2)] == [L(k,1), L(k,2)]
           k = k + 1;
        else
           POI(j,1) = Lc(i,1);
           POI(j,2) = Lc(i,2);
           j = j + 1;
        end
    end

    % add interface points to POI
    POI = [POI; I(:,1) I(:,2)];  
    POI = sortrows(POI,2);
    POI = sortrows(POI);

    % make Points_of_interest a matrix
    Points_of_interest = zeros(N,M);
    for k = 1:length(POI(:,1))
        Points_of_interest(POI(k,2),POI(k,1)) = 1;
    end
%}
    %{
    % check what points of interst look like
    surface = Points_of_interest;
    my_cmap = load(strcat(cp_pth,'aa_slainity.mat'));
    figure; pcolor(surface); hold on
    shading(gca,'interp')
    colormap(gray);
    %}

    %% STEP 4: Calculate the open angle for all POI

    thresholdimg = Water;
    precision = 360;
    [OA, OA_sparse] = Seaangles2(thresholdimg,precision); % Use the Seangles2 function provided by John Shaw 

    % find location within POI where the opening angle is 45 degrees 
    OA_map = full(OA_sparse);
    %POI_OA_map =  zeros([N,M]);
    shore_map =  zeros([N,M]);
    shore_points = [];

    %threshold = 135; %[DEG]
    for i = 1:N
        for j = 1:M
           if  OA_map(i,j) < threshold;
               shore_map(i,j) = 1;
               shore_points = [shore_points; X(i,j) Y(i,j)];
           end
        end
    end

    %Remove small connected components (e.g. ponds)
    CC = bwconncomp(shore_map==0); %takes all OA_map points outside set threshold (water points)
    shoreline_map_noponds = false(N,M);
    for i = 1:length(CC.PixelIdxList)
        num_pixels(i) = length(CC.PixelIdxList{i});
    end
    minlength = floor(max(num_pixels)/100);

    for i = 1:length(CC.PixelIdxList)
        if(length(CC.PixelIdxList{i})>minlength)
            shoreline_map_noponds(CC.PixelIdxList{i}) = true;
            i;
        end
    end
    mask = shoreline_map_noponds;
    shoreline_idx = (shoreline_map_noponds == 0); %switch back to land being true
    
    %elevation_mask = elevation;
    %elevation_mask(shoreline_map_noponds) = 1;
    
    %{
    % check what open angle map looks like without ponds
    figure;
    imagesc(shore_map); hold on
    imagesc(shoreline_map_noponds); hold on
    colormap(gray)
    %}

    [B] = bwboundaries(shoreline_idx);
    
    j = 1;
    shoreline_pixels = [];
    for k = 1:length(B)
        for i  = 1:length(B{k})
            if B{k}(i,1) ~= 1
                shoreline_pixels{k}(j,1) = B{k}(i,2); %B has swaped x and y axis
                shoreline_pixels{k}(j,2) = B{k}(i,1);
                j = j+1;
            end
        end
        % for plotting purposes, add buffer point before and after
        %shoreline_pixels{k} = [shoreline_pixels{k}(1,1), 1; shoreline_pixels{k}];
        %shoreline_pixels{k} = [shoreline_pixels{k};shoreline_pixels{k}(j,1), 1];
        j = 1;
    end
    
    %{
    figure
    imagesc(elevation); hold on
    scatter(shoreline_pixels{k}(:,1), shoreline_pixels{k}(:,2), 'o', 'r')
    %}

    
    
    %{
    for k = 1:length(B)
        boundary = [boundary; B{k}];
    end
    
    % remove pixels along the beach
    j = 1;
    shoreline_pixels = [];
    for i = 1:length(boundary(:,1))
        if boundary(i,1) ~= 1
            shoreline_pixels(j,1) = boundary(i,2);
            shoreline_pixels(j,2) = boundary(i,1);
            j = j+1;
        end
    end
    

    %}
    
end

