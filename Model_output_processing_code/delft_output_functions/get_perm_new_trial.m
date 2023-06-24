function [perm,hydcon,d10] = get_perm_new(input,sed_size,layer)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calculate permeability
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [N, M] = size(input,[1,2]);
    low_dmm = sed_size(end);
    high_dmm = sed_size(1);
    sed_num = length(sed_size);
    
    % Look at map of different sed types
    %{
    figure
    for i = 1:length(input(1,1,1,:))
        subplot(2,3,i)
        imagesc(input(:,:,1,i)*-1)
        colormap(gray(225));
        set(gca,'Xtick',[]); set(gca,'Ytick',[])
        set(gca,'YDir','normal')
    end
    %}
   
    
    % calculate total mass in each cell
    %msed_step = msed{step,1}; %gives 227 302 1 11 (M N layer sed) for a given time step
    totmass = zeros(N,M);
    for i = 1:N
        for j = 1:M
            totmass(i,j) = sum(input(i,j,1,:));
        end
    end 
    
    %{
    figure
    imagesc(totmass)
    set(gca,'YDir','normal')
    %}
    
    % calculate percent retained (also just the percent of mass for each sediment type)
    percent_mass = zeros([N,M,sed_num], 'double');
    for k = 1:sed_num
        for i = 1:N
            for j = 1:M
                %percent_mass(i,j,k) = input(i,j,layer,sed_num-k+1)/totmass(i,j); % [%] sed_num-k is to reverse (so it is small to large sediment)
                percent_mass(i,j,k) = input(i,j,layer,k)/totmass(i,j); % [%] sed_num-k is to reverse (so it is small to large sediment)
            end
        end
    end
    
    %{
    % Look at map of different sed types
    figure
    for i = 1:length(percent_mass(1,1,:))
        subplot(2,3,i)
        imagesc(percent_mass(:,:,i)*-1)
        colormap(gray(225));
        set(gca,'Xtick',[]); set(gca,'Ytick',[])
        set(gca,'YDir','normal')
    end
    %}

    % calculate cummulative percenter retained
    cumpercent_retained = zeros([N,M,sed_num], 'double');
    for i = 1:N
        for j = 1:M
            cumpercent_retained(i,j,:) = cumsum(percent_mass(i,j,:)); %[%]
        end
    end

    % calculate percent finer than 
    percent_finer = zeros([N,M,sed_num], 'double');
    for k = 1:sed_num
        for i = 1:N
            for j = 1:M
                percent_finer(i,j,k) = 1-cumpercent_retained(i,j,k); % [%]
            end
        end
    end

    % interpolate to find d50
    xnew = linspace(low_dmm, high_dmm, 100);  %sediment size d [mm] - more refined for interpolation
    
    d10 = nan(N,M);
    for i = 1:N
        for j = 1:M
            f2 = interp1(sed_size, reshape(percent_finer(i,j,:), [1,sed_num]), xnew, 'linear'); % linear interpolation
            [val,idx] = min(abs(f2-0.10)); %get index of value closest to 0.10
             if isnan(val)
                d10(i,j) = nan;
             else
                d10(i,j) = xnew(idx); % find x [d10 in mm] at y ~ 0.10
             end
             %{
             figure
             plot(xnew,f2,':.',xnew(idx),val,'x');
             set(gca, 'XDir','reverse')
            %}
        end
    end
    
    %{
    figure
    imagesc(d10)
    colorbar
    %set(gca,'Xtick',[]); set(gca,'Ytick',[])
    set(gca,'YDir','normal')
    %}
    
    % Hazen Methods
    C = 40;
    j = 2;
    d10_cm = d10/10;
    hazen = C.*(d10_cm.^j); %in cm/s
    hydcon = hazen/100; %in m/s 
    
    % Conversion factors: http://www.balleau.com/materials/EQUIVALENT_UNITS.PDF
    g = 9.81; %m/2
    rho = 1000; %in kg/m^3
    mu = 1.793E-3; % dynamic viscosity in kg/ms
    
    perm = (hydcon*mu)/(g*rho); %in m^s
    perm = perm/9.869233E-13 ; %in darcy: conversion from 1 m^2 = darcy
    
    %[m^2] conversion 1 darcy = 9.869233x10^-3 m^2
    %hydcon = perm * 5471277.189; % conversion from perm in m^2
    %perm = perm * 9.869233E-3; % converting perm to darcy
    
    %{
    figure
    imagesc(hydcon)
    colorbar
    set(gca,'ColorScale','log')
    set(gca,'YDir','normal')
    %}
    
    %{
    figure
    imagesc(perm)
    colorbar
    set(gca,'ColorScale','log')
    set(gca,'YDir','normal')
    %}
end

