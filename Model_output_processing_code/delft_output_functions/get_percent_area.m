function [AD_HP] = get_percent_area(perm_mask,AD,dx,low_threshold,high_threshold)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    highlow_perm = []; highPerm_area = 0;
    if exist('high_threshold','var') == 0
        for i = 1:length(perm_mask(:,1))
            for j = 1:length(perm_mask(1,:))
                if isnan(perm_mask(i,j)) == 1
                    highlow_perm(i,j) = 0;
                elseif perm_mask(i,j)  > low_threshold
                    highlow_perm(i,j) = 1;
                else
                    highlow_perm(i,j) = 0;
                end
            end
        end
    elseif exist('high_threshold','var') == 1
        for i = 1:length(perm_mask(:,1))
            for j = 1:length(perm_mask(1,:))
                if isnan(perm_mask(i,j)) == 1
                    highlow_perm(i,j) = 0;
                elseif perm_mask(i,j) > low_threshold && perm_mask(i,j) <= high_threshold
                    highlow_perm(i,j) = 1;
                else
                    highlow_perm(i,j) = 0;
                end
            end
        end
    end

    highPerm_area = get_area(highlow_perm,dx);
    AD_HP =  highPerm_area/AD;
end

