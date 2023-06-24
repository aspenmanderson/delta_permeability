function [perm,hydcon] = get_perm(input,sed_size,layer)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calculate permeability
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [N, M] = size(input,[1,2]);
    low_dmm = sed_size(1);
    high_dmm = sed_size(end);
    sed_num = length(sed_size);
    
    % calculate total mass in each cell
    %msed_step = msed{step,1}; %gives 227 302 1 11 (M N layer sed) for a given time step
    totmass = zeros(N,M);
    for i = 1:N
        for j = 1:M
            totmass(i,j) = sum(input(i,j,1,:));
        end
    end 

    % calculate percent retained (also just the percent of mass for each sediment type)
    percent_mass = zeros([N,M,sed_num], 'double');
    for k = 1:sed_num
        for i = 1:N
            for j = 1:M
                percent_mass(i,j,k) = input(i,j,layer,sed_num-k+1)/totmass(i,j); % [%] sed_num-k is to reverse (so it is small to large sediment)
            end
        end
    end

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
    % num_sed_type = sed_num;
    xnew = linspace(low_dmm, high_dmm, 100); %sediment size d [mm] - more refined for interpolation

    x50 = zeros(N,M);
    for i = 1:N
        for j = 1:M
            f2 = interp1(sed_size, reshape(percent_finer(i,j,:), [1,sed_num]), xnew, 'linear');  %'linear' 'cubic'
            [val,idx] = min(abs(f2-0.50)); %get index of value closest to 0.50
            x50(i,j) = xnew(idx); % find x [d50 in mm] at y ~ 0.50
        end
    end

    % Shepard Method
        %c and j values from Shepard for beach material
        %{
    c = 5.8087E-6; %12000 [gal_per_day/ft^2] /7.4805 * 3.6210E-9 -> 5.8087E-6 [cm^2] or 40.746 (48.89499801) for m
    j = 1.65; 
        %x50_m = x50/100; %[mm] -> [m]
    perm = c.*(x50.^j); % [m^2] create conductivity matrix for a single timestep 
        %saywer got values -3, -5, and -7 m/s 
    %perm_noout = rmoutliers(perm,'percentiles',[0.02 99.8]); %only used for scaling caxis
    %}
        
    c = 3500; %48.89499801; %12000 [gal_per_day/ft^2] * 0.003785 [gal to m^3] / 0.092903 [ft2 to m2] -> results in units of m
    j = 1.65; 
    x50_m = x50/100; %[mm] -> [m]
    shep = c.*(x50_m.^j); %  [gal_per_day/ft^2]
   
    % Conversion factors: http://www.balleau.com/materials/EQUIVALENT_UNITS.PDF
    g = 9.81; %m/2
    rho = 1000; %in kg/m^3
    mu = 1.793E-3; % dynamic viscosity in kg/ms
    
    hydcon = shep*4.7160E-7; %in m/2: conversion from 1 gpd/ft2 = 4.7160 x 10-5 cm/s = 4.7160 x 10-7 m/s
    perm = ((hydcon*mu)/(g*rho)); %in m^s
    perm = perm/9.869233E-13 ; %in darcy: conversion from 1 m^2 = darcy
    
    %[m^2] conversion 1 darcy = 9.869233x10^-3 m^2
    %hydcon = perm * 5471277.189; % conversion from perm in m^2
    %perm = perm * 9.869233E-3; % converting perm to darcy
    
    
end

