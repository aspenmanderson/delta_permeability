function [Water,L] = get_land(elev_clip, wl_clip)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    [N, M] = size(elev_clip); % resize grid without shore
    [X,Y] = meshgrid(1:1:M,1:1:N);  

    %P = [];
    %Land = zeros([N,M]);
    L = 0;
    Water = zeros([N,M]);
    %W = [];
    for i = 1:N
        for j = 1:M
            %P = [P; X(i,j) Y(i,j)];
            %if elev_clip(i,j) - wl_clip(i,j) < 0
            if elev_clip(i,j) < 0
               Water(i,j) = 1;
               %W = [W; X(i,j), Y(i,j)];
            %else
               %Land(i,j) = 1;
               %L = [L; X(i,j) Y(i,j)];
               %L = L + 1;
            end
        end
    end
    %{
    % need to remove land that builds up around the eastern and western boundaries (amonalies)
    CC = bwconncomp(Water==0); %takes all OA_map points outside set threshold (water points)
   
    R = [];
    for i = 1:length(CC.PixelIdxList)
        [R(:,1), R(:,2)] = ind2sub(size(Water), CC.PixelIdxList{i}); %give row and column of pixels
        %for j = 1:length(R)
            if sum(R(:,2) < 10) > 0 || sum(R(:,2) > M-10) > 0 %accounting for 5 pixels on both size of the boundary
                Water(CC.PixelIdxList{i}) = 1;
            end
        %end
        R = [];
    end
    %}
    % count land points
    L = sum(Water(:) == 0);
    
    
    
     %{
    % check what land looks like
    surface = Land;
   % my_cmap = load(strcat(cp_pth,'aa_slainity.mat'));
    figure; pcolor(surface); hold on
    shading(gca,'interp')
    colormap(gca,flipud(my_cmap.cm));
    %a = colorbar; ylabel(a,'Land','FontSize',18);

    % check what water looks like
    surface = Water;
   % my_cmap = load(strcat(cp_pth,'aa_slainity.mat'));
    figure; pcolor(surface);
    shading(gca,'interp')
    colormap(gca,flipud(my_cmap.cm));
    %a = colorbar; ylabel(a,'Water density (kg/m^3)','FontSize',18);
    %}
end

