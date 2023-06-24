function [clipped_map, N_clipped, M_clipped] = clip_map(map, clip_dim)

% Clips a map to a specified size and returns the new dimensions
    % For Delft3D outputs, clipping the domain:
         % - removes sediment accumulation along the boundary
         % - removes the shoreline

[N, M] = size(map); 
clipped_map = map(clip_dim(3):N-clip_dim(1),clip_dim(4):M-clip_dim(2));
[N_clipped, M_clipped] = size(clipped_map); 

end

