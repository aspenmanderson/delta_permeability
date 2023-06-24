function [sinuosity] = get_shoreline_rugosity(shoreline_pixels,dx)
%UNTITLED2 Summary of this function goes here
%   Sinuosity is calculated by dividing the length of the shore by the 
%   straight line distance between the end points

    x  = [];
    y_rough  = [];
    for k = 1:length(shoreline_pixels)
        x = [x; shoreline_pixels{k}(:,1)];
        y_rough = [y_rough; shoreline_pixels{k}(:,2)];
    end

    y_smooth = smoothdata(y_rough,'gaussian',30); %smoothing window = 30

    %{
    % check what shorelines look like
    figure
    plot(x,y_rough); hold on
    plot(x,y_smooth)
    legend('rough shoreline','smoothed shoreline')
    %}

    % calculate the distance between each point along the line
    dy = dx; %assuming square grid cell
    for i = 1:length(y_rough)-1
        distance_rough(i) = sqrt((x(i+1)-x(i)).^2 + (y_rough(i+1)-y_rough(i)).^2)*sqrt(dx^2 + dy^2);
            % sqrt(dx^2 + dy^2) assumes distances are meassured in the middle
            % of each cell
    end

    tot_dist_rough = sum(distance_rough);

    for i = 1:length(y_smooth)-1
        distance_smooth(i) = sqrt((x(i+1)-x(i)).^2 + (y_smooth(i+1)-y_smooth(i)).^2)*sqrt(dx^2 + dy^2);
    end

    tot_dist_smooth = sum(distance_smooth);


    % calculate sinuosity
    sinuosity = tot_dist_rough/tot_dist_smooth;

end

