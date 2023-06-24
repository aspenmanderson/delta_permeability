function [shape, len] = get_shape(elevation)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    [N, M] = size(elevation); % resize grid without shore

    % find maximum length (N direction)
    cw = [];
    for i = 1:N
        cw(i) = nnz(elevation(i,:));
    end

    wid = max(cw);

    % find maximum length (M direction)
    cl = [];
    for j = 1:M
        cl(j) = nnz(elevation(:,j));
    end

    len = max(cl);

    % calculate shape
    shape = wid/(2*len);

end

