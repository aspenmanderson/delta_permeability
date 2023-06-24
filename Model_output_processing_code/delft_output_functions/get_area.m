function [area] = get_area(input, dx)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

[N, M] = size(input); % resize grid without shore


k = 0;
for i = 1:N
    for j = 1:M
        if input(i,j) ~= 0
            k = k + 1;
        end
    end
end

dy = dx; % assuming square pixels
area_m2 = k*dx*dy; % number of land pixels times the size of the pixel [M^2]
area = area_m2*1e-6; % converting to [km^2]

end

