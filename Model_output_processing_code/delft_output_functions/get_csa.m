function [csa] = get_csa(elevation, thickness)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    % find width along initial shoreline
    wid = nnz(elevation(1,:));

    % find mean thickness along initial shoreline
    depth = mean(thickness(1,:));

    csa = wid*depth;
    
end

