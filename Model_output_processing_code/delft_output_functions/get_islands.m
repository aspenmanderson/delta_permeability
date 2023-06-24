function [num_islands,num_plain] = get_islands(mask,AD,dx)
%count number of islandse
%   Detailed explanation goes here

    CC = bwconncomp(mask);
    num_objects  = CC.NumObjects;
    %{
    num_islands = 0;
    for i = 1:length(CC.ImageSize)
        if CC.ImageSize(i) < AD/dx/dx/10 % islands have to be less than 10 percent total delta area 
            num_islands = num_islands + 1;
        end
    end
    
    num_plain = num_objects - num_islands;
    %}
    if num_objects == 1
        num_plain = num_objects;
        num_islands = 0;
    else
        num_plain = 1;
        num_islands = num_plain - num_objects;
    end
    
    
end

