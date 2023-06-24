function [Dgrd] = get_gradient(shoreline_pixels, water_level, dx )
%Computes the average topset gradient
%   we measure rays from the delta apex to the shoreline. The elevation
%   of the delta apex is defined by averaging all water surface elevations 
%   within a 500m swath around the head of the delta.

    % find the middle of the main channel (middle of domain 
    [N, M] = size(water_level); % resize grid without shore
    [X,Y] = meshgrid(1:1:M,1:1:N);
    center_x = ceil(N/2);
    center_y = 1;
    
    grad = [];
    for k = 1:length(shoreline_pixels)
        for i = 1:length(shoreline_pixels{k}(:,1))
            diff(i) = water_level(center_x,center_y) - water_level(shoreline_pixels{k}(i,2),shoreline_pixels{k}(i,1));
            distance(i) = sqrt((center_x-shoreline_pixels{k}(i,2))^2 + (center_y-shoreline_pixels{k}(i,1)*dx)^2);
        end
        grad = [grad; diff(:)./distance(:)];
    end
    
    Dgrd = mean(grad)*-1; % swap sign so a gradient toward the shoreline is positive
    
end