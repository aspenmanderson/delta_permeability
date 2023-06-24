function [masked_map] = mask_map(masked_map,mask)

% Returns a map with only piixels where mask = true
    % returned map is the same size as the original map
    
masked_map(mask) = nan;

%{
figure
imasgesc(masked_map)
colorbar
%}

end

