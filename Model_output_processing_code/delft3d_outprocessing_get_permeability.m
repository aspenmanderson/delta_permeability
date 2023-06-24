% Author: Aspen Anderson
% This code processes Delft3D outputs using a variety of developed
% functions
% Last Modified: 01/06/2023

clc; close all; clear all;

% add function files
addpath("C:\Delft3d\repo\src\tools_lgpl\matlab\quickplot\progsrc"); %D3D
addpath("C:\Users\labuser10\OneDrive - sfu.ca\Working\SFU\Research\Modeling\delta_modeling_paper_codes\output_processing\delft_output_functions"); %functions to process outputs

%cp_pth = ("C:\Program Files\Deltares\Delft3D 4.03.01\win64\quickplot\bin\colormaps\"); %colorplot location

set(0,'DefaultAxesFontSize', 18)
set(0,'DefaultAxesFontName','Times') 
set(0,'DefaultAxesFontWeight','Demi')
set(0,'DefaultAxesLineWidth',2)

%%
tic
for kk = 1
    folder = strcat('I:\Aspen_Delft3D_local\Summer_2020_trials\output_CC_try6_REDO\');
    ex_pth = folder;
    
    % need to generate all Delft3D output before running. This grabs all output files (.dat)
    if kk == 1
        model_type = 'F';
        ran = 1:75;
    elseif kk == 2
        model_type = 'S';
        ran = [1:39 41:51 53:56];
    elseif kk == 3
        model_type = 'SS';
        ran = [1:17 19:20 22:34 36:41 43:45 47 48 50:52 54:62 64:77 79:87 89 90 92:94 96 97];
    end
    
    for mod_num = ran
        
        close all
        
        % define model input
        model = strcat(model_type,num2str(mod_num))
        input = strcat(folder, model,'.dat'); %trim-trial
        
        % resetting end_step for models that interfere with boundaries
        clear end_step
        if kk == 3 % only SS models 
            if mod_num == 76 || mod_num == 96 || mod_num == 97
                end_step = 10;
            elseif mod_num == 20
                end_step = 20;
            elseif mod_num == 77 || mod_num == 89
                end_step = 30;
            elseif mod_num == 19 || mod_num == 27 
                end_step = 40;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Load data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % open .dat output file from D3D
        nfs = vs_use(input); % opens nefis file (d3d output)
        vs_disp(nfs) % list of groups and varibales contianed in nefis files
            %how to access Data2{24, 1}

        % get grid params 
        dx = 200; % size of grid cell in x direction [M]
        dy = 200; % size of grid cell in y direction [N]

        M = vs_get(nfs, 'map-const', '', 'MMAX', ''); % number of m-grid points
        N = vs_get(nfs, 'map-const', '', 'NMAX', ''); % number of n-grid points
        K = vs_get(nfs, 'map-const', '', 'KMAX', ''); % number of layers points

        % get time params
        times = vs_get(nfs, 'map-const', '', 'TUNIT', ''); %  time related to seconds [S]
        itavgs = vs_get(nfs, 'map-infavg-serie', '', 'ITAVGS', '');  %end time [S]  
        itdate = vs_get(nfs, 'map-const', '', 'ITDATE', ''); % initial date and time [YYYYMMDD]
        dt = vs_get(nfs, 'map-const', '', 'DT', ''); % timestep [S]
        morft = vs_get(nfs, 'map-infsed-serie', '', 'MORFT', ''); % morphological time [DAYS since simulation started]
        morfac = vs_get(nfs, 'map-infsed-serie', '', 'MORFAC', ''); % morphological time [DAYS since simulation started]
        morfac = morfac{1};

        % get flow params
        waterlevel = vs_get(nfs, 'map-series', '', 'S1', ''); % water-level in zeta point [M]

        % get sed params
        sed_num = vs_get(nfs, 'map-const', '', 'LSTCI', ''); % number of layers points
        sed_num = sed_num - 1; %takes out the value that is just the total/cummulative
        msed = vs_get(nfs, 'map-sed-series', '', 'MSED', ''); % mass of sediment in layer [KG/M3]
        lyrfrac = vs_get(nfs, 'map-sed-series', '', 'LYRFRAC', ''); % volume fraction of sediment in layer [-]
        dp_bedlyr = vs_get(nfs, 'map-sed-series', '', 'DP_BEDLYR', ''); % vertical position of sediment layer interface [M]
        namsed =  vs_get(nfs, 'map-const', '', 'NAMSED', ''); % name of sediment fraction
        thick  = vs_get(nfs, 'map-const', '', 'THICK', ''); % fraction part of layer thickness of total water-height
        dps = vs_get(nfs, 'map-sed-series', '', 'DPS', ''); % bottom depth (zeta point)
        DP0 = vs_get(nfs, 'map-const', '', 'DP0', ''); % initital bottom depth
        rca = vs_get(nfs, 'map-sed-series', '', 'RCA', ''); % mass of sediment in layer [KG/M3]

        namcon = vs_get(nfs, 'map-const', '', 'NAMCON', ''); % name of concentrations
        r1 = vs_get(nfs, 'map-series', '', 'R1', ''); % index 1 = salinity concentrations
        
        sed_size = [2, .640, .125, 0.0313, 0.0078, 0.0020];
        layer = 1;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Initialize parameters and calculate values for each timestep
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        AD = [];  CN = []; Tcl = []; L = []; Dsh = [];
        PERM = []; AD_HP = []; CSA = []; q = []; q_min = []; q_max = []; HYDCON = []; 
        sill = []; range = [];
        clip_dim = [200, 3, 29, 3]; % [north, east, south, west]             
                % removes the beach (54 pixels along southern boundary) and boundary influences (26 pixels along north and 3 around east & west boundaries)

        start = 10 ; % skip spin up time
        if exist('end_step','var') == 0 % if variable is not previously defined
            end_step = length(waterlevel); 
        end
        output_times = start:3:end_step; 


        ii = 1;
        for step = output_times(end)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Format outputs
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % get elevation  
            elevation = (dps{step, 1} - waterlevel{step, 1})*-1;         
            wl = waterlevel{step, 1};
            [N, M] = size(elevation);


            [elev_clip, N_clip, M_clip] = clip_map(elevation, clip_dim);  
            [DP0_clip] = clip_map(DP0, clip_dim);
            [wl_clip] = clip_map(wl, clip_dim);
            wh_clip = wl_clip-elev_clip;
            [bd_clip] = clip_map(DP0,clip_dim);
            thick = elev_clip - bd_clip*-1;

            % get mask
                % distinguests land, water, and shoreline pixels
                % mask == True where land pixels are present 
            [Water,L_points] = get_land(elev_clip, wl_clip);    % check to see if there is land yet

            if L_points > 375 % enough pixels to calculate a convex hull (stating condition- equals 5 km^2, same as smallest deltas from Lit)
                %condition was 10 before change on Jan 18, 2021, worked well
                %at some point I calculated 125 cells is 5 km^2 (375 is 15
                %km^2)

                [shore_mask, shoreline_pixels] = get_shoreline(Water,45); %threshold = 135 [DEG]. 45 works better on the deltas that come close to the western boundary

                % apply land_mask to elevation map
                [elev_mask] = mask_map(elev_clip,shore_mask);
                %[DP0_mask] = mask_map(DP0_clip,shore_mask);
                [wl_mask] = mask_map(wl_clip,shore_mask);
                %[wh_mask] = mask_map(wh_clip,shore_mask);



                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Calculate variables
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % area
                [AD(ii)] = get_area(elev_mask,dx);

                % get permeability map and format
                [perm,hydcon,d10] = get_perm_hazen(msed{step,1},sed_size,layer);
                %[perm,hydcon,d10] = get_perm_USBR(msed{step,1},sed_size,layer);
                [perm_clip] = clip_map(perm, clip_dim); % clip domain and apply mask to isolate delta
                [perm_mask] = mask_map(perm_clip,shore_mask);
                [hydcon_clip] = clip_map(hydcon, clip_dim); % clip domain and apply mask to isolate delta
                [hydcon_mask] = mask_map(hydcon_clip,shore_mask);

                % remove outliers within the delta
                   % - need to remove 0's otherwise skews the geomeans
                perm_noout = filloutliers(perm_mask,'clip','percentiles',[0.01 99.9]);
                [perm_mask] = mask_map_nan(perm_noout,shore_mask); 

                hydcon_noout = filloutliers(hydcon_mask,'clip','percentiles',[0.01 99.9]);
                [hydcon_mask] = mask_map_nan(hydcon_noout,shore_mask); 
                
                [d10_mask] = mask_map_nan(d10,shore_mask); 

                % geomean of permeability 
                [PERM(ii)] = geomean(reshape(perm_mask,[1,size(perm_mask,1)*size(perm_mask,2)]),'omitnan'); % convert matrix to vector find geomean
                [HYDCON(ii)] = geomean(reshape(hydcon_mask,[1,size(hydcon_mask,1)*size(hydcon_mask,2)]),'omitnan'); % convert matrix to vector find geomean

                % percent of high permeability area
                [AD_HP(ii)] = get_percent_area(perm_mask,AD(ii),dx,30); %geomean of the delta
                [AD_HP_csand(ii)] = get_percent_area(hydcon_mask,AD(ii),dx,5E-5,6E-3); 
                [AD_HP_fsand(ii)] = get_percent_area(hydcon_mask,AD(ii),dx,2E-6,5E-5); 
                [AD_HP_silt(ii)] = get_percent_area(hydcon_mask,AD(ii),dx,1E-9,2E-6);
                    % values from Domenico and Schwarts, 1990: http://www.aqtesolv.com/aquifer-tests/aquifer_properties.htm

                % channel network properties
                [shape(ii),len] = get_shape(elev_mask); 
                [CN(ii), Tcl(ii), L(ii), Ncw(ii), bw_centerline_noponds, channel_width_map] = get_rivers_new(elev_mask*-1,shoreline_pixels,dx,len);

                % Comparison of channel network and high perm areas 
                % connectivity of high permeable areas
                
                thrshold(ii) = prctile(reshape(perm_mask,[1,size(perm_mask,1)*size(perm_mask,2)]),75); %90th percentile of the perm in that delta
                [C(ii),NumC(ii),SizeC(ii),percent_hp_and_cn(ii)] = get_connectivity(perm_mask,thrshold(ii),bw_centerline_noponds, AD(ii), dx, elev_clip); %90th percentile of the geomeans
                
                % delta plain properties
                [SLOPE(ii),grad] = get_gradient2(shoreline_pixels, elev_mask, dx);
                [Dsh(ii),clino_profile,x_clino] = get_subaq_clino2(shore_mask, elev_clip, DP0_clip);
                [shore_rug(ii)] = get_shoreline_rugosity(shoreline_pixels,dx);
                
                    % A=1 are semicircular
                    % A<1 are elongate with long axes perpendicular to the beach
                    % A>1 have long axes that are beach parallel

                %[num_islands(ii),num_plain(ii)] = get_islands(shore_mask,AD,dx);
                    % num plain should always be 1 if shoreline is calculated correctly

                %cf = 5.3997e+04; % [m] conversion factor: 5471277.189*9.869233E-3;
                q(ii) = HYDCON(ii)*SLOPE(ii); % [m3/s]
                q_min(ii) =  min(min(hydcon_mask))*min(abs(grad))*-1;
                q_max(ii) =  max(max(hydcon_mask))*max(abs(grad))*-1;


                time_step(ii) = step;
                morpho_time(ii) = morft{step}/360;


                ii = ii + 1; 
            end
        end

       % Export figure for last time step
       % export_delft_figure(M_clip,N_clip,elev_clip,morft{(step)},ex_pth,model, 0)  %last argument: 0 = no label, 1 = label
        if L_points > 375 
           % export_perm_figure(M_clip,N_clip,perm_mask,morft{(step)},ex_pth,strcat(model,'_perm'), 0) 
            %scatter(shoreline_pixels{K}(:,1), shoreline_pixels{K}(:,2), '.', 'r')
           % export_perm_figure(M_clip,N_clip,hydcon_clip,morft{(step)},ex_pth,strcat(model,'_hydcon'), 0)
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Export variables to Excel
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        if sum(CN) > 0 && L_points > 375 && shape(end) < 5 % shape < 5 to make sure shoreline protrudes (length has to be greater than 1/10th of the width)
            % spread sheet for each delta model
            %T = table(time_step', morpho_time', AD', AD_high', PERM', CN', Tcl', L', Dsh', shape', shore_rug', SLOPE');
            %T.Properties.VariableNames = {'Time step', 'Morphological time', 'AD', 'AD_high', 'Perm', 'CN', 'Tcl', 'L', 'Dsh', 'shape', 'Shoreline rugoisty','SLOPE'};
            %writetable(T,strcat(ex_pth,model,"output_statistsics.xlsx"));

            % spread sheet for all delta models at last time step
            data_table = [mod_num time_step(end) morpho_time(end) AD(end) CN(end) Tcl(end)  L(end) Ncw(end) Dsh(end) shape(end) shore_rug(end) PERM(end) HYDCON(end) geomean(geomean(d10_mask,'omitnan'),'omitnan') C(end) NumC(end) SizeC(end) thrshold(end) percent_hp_and_cn(end) SLOPE(end) q(end) q_min(end) q_max(end) AD_HP(end) AD_HP_csand(end) AD_HP_fsand(end) AD_HP_silt(end) range(1:end) sill(1:end)]; 
            table_name = strcat(ex_pth,"ALL_",model_type,"models_output_statistics.xlsx");
            xlswrite(table_name,data_table,'Sheet1',['A' num2str(mod_num+1)]);
        end
    end
    
end
toc
