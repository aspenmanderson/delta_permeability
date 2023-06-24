function [Water,W] = get_water(elev_clip, wl_clip)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    [N, M] = size(elev_clip); % resize grid without shore
    [X,Y] = meshgrid(1:1:M,1:1:N);  

    %P = [];
    %Land = zeros([N,M]);
    %L = [];
    Water = zeros([N,M]);
    W = [];
    for i = 1:N
        for j = 1:M
            %P = [P; X(i,j) Y(i,j)];
            if elev_clip(i,j) - wl_clip(i,j) > 0
               %Land(i,j) = 1;
               %L = [L; X(i,j) Y(i,j)];
            else
               Water(i,j) = 1;
               W = [W; X(i,j), Y(i,j)];
            end
        end
    end
    
     %{
    % check what land looks like
    surface = Land;
    my_cmap = load(strcat(cp_pth,'aa_slainity.mat'));
    figure; pcolor(surface); hold on
    shading(gca,'interp')
    colormap(gca,flipud(my_cmap.cm));
    %a = colorbar; ylabel(a,'Land','FontSize',18);

    % check what water looks like
    surface = Water;
    my_cmap = load(strcat(cp_pth,'aa_slainity.mat'));
    figure; pcolor(surface);
    shading(gca,'interp')
    colormap(gca,flipud(my_cmap.cm));
    %a = colorbar; ylabel(a,'Water density (kg/m^3)','FontSize',18);
    %}
end

