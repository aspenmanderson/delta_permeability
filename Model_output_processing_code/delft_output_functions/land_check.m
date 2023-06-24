function [P, Land, L, Water] = land_check(elev_clip, wl_clip)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    [N, M] = size(elev_clip); % resize grid without shore
    [X,Y] = meshgrid(1:1:M,1:1:N);  

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
end

