function [Dgrd,grad] = get_gradient2(shoreline_pixels, elevation, dx )
%Computes the average topset gradient
%   we measure rays from the delta apex to the shoreline. The elevation
%   of the delta apex is defined by averaging all water surface elevations 
%   within a 500m swath around the head of the delta.

% using water level mkes gradient go up and down for tidal deltas. Need to
% use sediment thinkness

    % find the middle of the main channel (middle of domain 
    [N, M] = size(elevation); % resize grid without shore
    [X,Y] = meshgrid(1:1:M,1:1:N);
    
    % for x center, using max elevation within the 10 cells on each side of the channel
        % incase there is erosion along the initial shore due to channel formation, as observed in some
        % models
    center_y = 1;
    % center_x1 = 195; % the two pixels on the side of the channel
    % center_x2 = 207;
    
    grad = [];
    for k = 1:length(shoreline_pixels)
        for i = 1:length(shoreline_pixels{k}(:,1))
            %diff(i) = ((elevation(center_y,center_x1)+elevation(center_y,center_x2))/2) - elevation(shoreline_pixels{k}(i,2),shoreline_pixels{k}(i,1));
            %distance(i) = sqrt((((elevation(center_y,center_x1)+elevation(center_y,center_x2))/2)-shoreline_pixels{k}(i,2))^2 + (center_y-shoreline_pixels{k}(i,1)*dx)^2);
            diff(i) = ((max(elevation(center_y,185:195))+max(elevation(center_y,207:217)))/2) - elevation(shoreline_pixels{k}(i,2),shoreline_pixels{k}(i,1));
            distance(i) = sqrt((((max(elevation(center_y,207:217))+max(elevation(center_y,207:217)))/2)-shoreline_pixels{k}(i,2))^2 + (center_y-shoreline_pixels{k}(i,1)*dx)^2);
        end
        grad = [grad; diff(:)./distance(:)]*-1;
    end
    
    Dgrd = mean(grad); % swap sign so a gradient toward the shoreline is positive
    
end