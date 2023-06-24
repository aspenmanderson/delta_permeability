%% Distance-Based Generalized Sensitivity Analysis
% Author: Celine Scheidt
% Date: August 2012
% Updated: January 2014

clc; clear all; close all;

% add function files
cp_pth = ("C:\Program Files\Deltares\Delft3D 4.03.01\win64\quickplot\bin\colormaps\"); %colorplot location

set(0,'DefaultAxesFontSize', 18)
set(0,'DefaultAxesFontName','Times') 
set(0,'DefaultAxesFontWeight','Demi')
set(0,'DefaultAxesLineWidth',2)

% Set parameters for run
threshold_f = 0.5; threshold_w = 0.4; threshold_t = 0.35;% for normalized force ratio threshold
range = '2:229'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

in_pth = ("C:\Users\labuser10\OneDrive - sfu.ca\Working\SFU\Research\Modeling\delta_modeling_paper_codes\output_processing\"); 
data_name = ("MASTER_all_inputs_and_outputs.xlsx");

% import data
opts = detectImportOptions(strcat(in_pth,data_name));
opts.Sheet = 'data';

% select delta name
opts.SelectedVariableNames = [1]; 
opts.DataRange = range; %Setting range of values 2:75 = F, 76:101 = S, 102:215 = SS 
[model_name] = readvars(strcat(in_pth,data_name),opts);
model_num = 1:1:length(model_name); 
N = length(model_num);

% select input parameters
opts.SelectedVariableNames = [4 5 6 7 8 9 10 11 12 13 14]; 
[discharge, wave, tide, nQr, nQw, nQt, sed_con, bath_grad, river_width, sed_dis, dmm] = readvars(strcat(in_pth,data_name),opts);

opts.SelectedVariableNames = [17 18 20 22 23]; 
[cn, tcl, ncw, shape, shore_rug] = readvars(strcat(in_pth,data_name),opts);

% select output parameters
opts.SelectedVariableNames = [3 16 24 25 26 32 33]; % 16 17 21 24]; 
[SedFormed, AD, PERM, HYDCON, d10, SLOPE, q] = readvars(strcat(in_pth,data_name),opts);
SLOPE = abs(SLOPE); %change for fit distribution later in code
q = abs(q);

geomean(PERM,'omitnan')

% select output parameters
opts.SelectedVariableNames = [27 28 29 30 31 37 38 39]; % 16 17 21 24]; 
[C, NumC, SizeC, threshold, percent_hp_cn, AD_HP_high, AD_HP_med, AD_HP_low] = readvars(strcat(in_pth,data_name),opts);
SLOPE = abs(SLOPE); %change for fit distribution later in code

% range variables
opts.SelectedVariableNames = [40:45]; 
[a b c d e f] = readvars(strcat(in_pth,data_name),opts);
range = [a, b, c, d, e, f];

% sill variables
opts.SelectedVariableNames = [46:51]; 
[a b c d e f] = readvars(strcat(in_pth,data_name),opts);
sill = [a, b, c, d, e, f];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Categorize models by force distribution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:N
    if floor(discharge(i)) == 1845
        Pf(i,1) = 0.5;
    elseif floor(discharge(i)) == 6768
        Pf(i,1) = 0.75;
    else
        Pf(i,1) = 0.9;
    end
    
    if wave(i) == 0
        Pw(i,1) = 0;
    elseif round(wave(i),1) == 1.0
        Pw(i,1) = 0.25;
    elseif round(wave(i),1) == 1.4
        Pw(i,1) = 0.50;
    elseif round(wave(i),1) == 1.9
        Pw(i,1) = 0.75;
    else
        Pw(i,1) = 0.9;
    end
    
    if tide(i) == 0
        Pt(i,1) = 0;
    elseif round(tide(i),1) == 0.4
        Pt(i,1) = 0.25;
    elseif round(tide(i),1) == 1.1
        Pt(i,1) = 0.50;
    elseif round(tide(i),1) == 2.1
        Pt(i,1) = 0.75;
    else
        Pt(i,1) = 0.9;
    end 
end

nFRf = Pf./(Pf+Pw+Pt);
nFRw = Pw./(Pf+Pw+Pt);
nFRt = Pt./(Pf+Pw+Pt);
nFR_sum = nFRf+nFRw+nFRt; %should be 1 for all cases

Pw_Pr = Pw./Pf;
Pt_Pr = Pt./Pf;
Pw_Pt = Pw./Pt;

for i = 1:N
    if isnan(Pw_Pt(i,1)) == 1 || Pw_Pt(i,1) == Inf
        Pw_Pt(i,1) = 0;
    end
end

TF = (Pf+Pw+Pt)/3; %Fluxuates from 0 to 3

idx_all = zeros(N,1); idx_fluvial = zeros(N,1); idx_wave = zeros(N,1); idx_tidal = zeros(N,1); idx_mixed = zeros(N,1);
for i = 1:N
    if SedFormed(i,1) == 1 && isnan(AD(i,1)) == 0 && isnan(SedFormed(i,1)) == 0
        idx_all(i,1) = 1;
        if nFRf(i,1) >= threshold_f 
            idx_fluvial(i,1) = 1;
        end
        if nFRw(i,1) >= threshold_w
            idx_wave(i,1) = 1;
        end
        if nFRt(i,1) >= threshold_t
            idx_tidal(i,1) = 1; 
        end
        if nFRf(i,1) <= threshold_f &&  nFRw(i,1) <= threshold_w && nFRt(i,1) <= threshold_t
            idx_mixed(i,1) = 1;
        end
    end
end

length(idx_fluvial(idx_fluvial == 1))
length(idx_wave(idx_wave == 1))
length(idx_tidal(idx_tidal == 1))
length(idx_mixed(idx_mixed == 1))
length(idx_all(idx_all == 1))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Groundwater parameters: distrubitions and export table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculating variability in mean q for different delta types
for ii = [5 1:3]
    
    % select model responce data
    if ii == 1
        idx = idx_fluvial;
        length(idx_fluvial(idx_fluvial == 1))
    elseif ii == 2
        idx = idx_wave;
        length(idx_wave(idx_wave == 1))
    elseif ii == 3
        idx = idx_tidal;
        length(idx_tidal(idx_tidal == 1))
    else
        idx = idx_all; % all models
    end

end

ex_pth = "C:\Users\labuser10\OneDrive - sfu.ca\Working\SFU\Research\Modeling\conductivity_paper\dgsa_perm\hydro_parameter_ranges"; %output destination


input_idx = [idx_fluvial, idx_wave, idx_tidal, idx_all];
input_var = [PERM, HYDCON,  SLOPE, q, AD_HP_high, AD_HP_med, AD_HP_low];
input_var_name = {'PERM'; 'HYDCON'; 'SLOPE'; 'q'; 'AD_HP_sand'; 'AD_HP_silt'; 'AD_HP_clay'};
input_unit_name = {''; ''; ''; ''; ''; ''};
find_distribution(input_var,input_idx,input_var_name,input_unit_name,ex_pth,1)

input_var = [AD_HP_med];
input_var_name = {'AD_HP_silt'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make pie charts for median percent area within each type of delta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure 

input_idx = []; input_var = [];
input_idx = [idx_fluvial, idx_wave, idx_tidal];
input_var = [AD_HP_high, AD_HP_med, AD_HP_low];
labels = {'Sand'; 'Silt'; 'Clay'};

subplot(1,4,1)
pie([mean(AD_HP_high(idx_fluvial == 1)) mean(AD_HP_med(idx_fluvial == 1)) mean(AD_HP_low(idx_fluvial == 1))]);
subplot(1,4,2)
pie([mean(AD_HP_high(idx_wave == 1)) mean(AD_HP_med(idx_wave == 1)) mean(AD_HP_low(idx_wave == 1))]);
subplot(1,4,3)
pie([mean(AD_HP_high(idx_tidal == 1)) mean(AD_HP_med(idx_tidal == 1)) mean(AD_HP_low(idx_tidal == 1))]);
subplot(1,4,4)
pie([nanmean(AD_HP_high) nanmean(AD_HP_med) nanmean(AD_HP_low)]);
legend('high','med','low');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make box plots for CN, Tcl, Shape, Sr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

input_idx = []; input_var = [];
input_idx = [idx_fluvial, idx_wave, idx_tidal];
input_var = [cn, shape, shore_rug];
label = {'Number of Channels'; 'Shape'; 'Shoreline Rugosity'};


for i = 1:length(input_var(1,:))
    figure
    %input_var = [];

    % data in cell array
    datac{1} = input_var(idx_fluvial == 1,:); % fluvial
    datac{2} = input_var(idx_wave == 1,:); % wave
    datac{3} = input_var(idx_tidal == 1,:);  % tidal
    datac{4} = input_var(:,:);
    
    group_names = {'Fluvial';'Wave';'Tidal';'All'};
    condition_names = {'0';'30';'60';'90';'120';'150'};

    % an alternative color scheme for some plots
    c =  [0.45, 0.80, 0.69;...
          0.98, 0.40, 0.35;...
          0.55, 0.60, 0.79;...
          0.90, 0.70, 0.30]; 

    h = daboxplot(datac,'xtlabels',condition_names,'colors',c,'whiskers',0,'scatter',1,'scattersize',15,'scatteralpha',0.5,'boxspacing',0.8); 
   % ylabel('Performance');
     % set x axis scale
     %xlim([min(input_var(:,i)) max(input_var(:,i))]);
    %if i == 2 
    %    set(gca,'YScale','log');
    %    ylim([1E-1 2E4]);
    %else
    %    ylim([0 2.3E4]);
    %end
    ylabel(label{i})
    set(gca,'FontSize',9.5);
   % xl = xlim; 
   % xlim([xl(1), xl(2)+0.2]);    % make more space for the legend
   % legend([h.bx(1,:)],group_names);

end
                                                                     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make sticker plots for CN, Tcl, Shape, Sr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

input_idx = []; input_var = [];
input_idx = [idx_fluvial, idx_wave, idx_tidal];
input_var = [cn, shape, shore_rug];
label = {'Number of Channels'; 'Shape'; 'Shoreline Rugosity'};



for i = 1:length(input_var(1,:))   
    %boxplot([input_var(idx_fluvial == 1,i), input_var(idx_wave == 1,i), input_var(idx_tidal == 1,i), input_var],)

    figure
    % data in a numreic array (+ grouping indices)
    data = NaN(length(input_var(:,i)),4);
    data(1:length(input_var(idx_fluvial == 1,i)),1) = input_var(idx_fluvial == 1,i);
    data(1:length(input_var(idx_wave == 1,i)),2) = input_var(idx_wave == 1,i);
    data(1:length(input_var(idx_tidal == 1,i)),3) = input_var(idx_tidal == 1,i);
    data(:,4) = input_var(:,i);
    %group_inx = [ones(1,30), 2.*ones(1,30), 3.*ones(1,30), 4.*ones(1,30)];

    %group_names = {'Humans', 'Dogs' , 'God', 'Potato'};
    condition_names = {'Fluvial';'Wave';'Tidal';'All'};

    % an alternative color scheme for some plots
    c =  [0.45, 0.80, 0.69;...
          0.98, 0.40, 0.35;...
          0.55, 0.60, 0.79;...
          0.90, 0.70, 0.30]; 

    h = daboxplot(data,'xtlabels',condition_names,'colors',c,'whiskers',0,'scatter',1,'scattersize',25,'scatteralpha',0.5,'boxspacing',0.2); 
   % ylabel('Performance');
     % set x axis scale
    ylabel(label{i})
    set(gca,'FontSize',9.5);
    xl = xlim; 

end

ex_pth = "C:\Users\labuser10\OneDrive - sfu.ca\Working\SFU\Research\Modeling\conductivity_paper\dgsa_perm\delta_parameter_ranges"; %output destination

input_var_name = {'Number of Channels'; 'Shape'; 'Shoreline Rugosity'};
input_unit_name = {''; ''; '';};
find_distribution(input_var,input_idx,input_var_name,input_unit_name,ex_pth,1)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make sticker plots for conductivity, number of conductive areas, size of connective body and percent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

input_idx = []; input_var = [];
input_idx = [idx_fluvial, idx_wave, idx_tidal];
input_var = [C, NumC, SizeC, percent_hp_cn];
label = {'Connectivity'; 'Number of Connective Bodies'; 'Size of Largest Connective Bodies'; 'Percent of Cn that is HP'};

for i = 1:length(input_var(1,:))   
    %boxplot([input_var(idx_fluvial == 1,i), input_var(idx_wave == 1,i), input_var(idx_tidal == 1,i), input_var],)

    figure
    % data in a numreic array (+ grouping indices)
    data = NaN(length(input_var(:,i)),4);
    data(1:length(input_var(idx_fluvial == 1,i)),1) = input_var(idx_fluvial == 1,i);
    data(1:length(input_var(idx_wave == 1,i)),2) = input_var(idx_wave == 1,i);
    data(1:length(input_var(idx_tidal == 1,i)),3) = input_var(idx_tidal == 1,i);
    data(:,4) = input_var(:,i);
    %group_inx = [ones(1,30), 2.*ones(1,30), 3.*ones(1,30), 4.*ones(1,30)];

    %group_names = {'Humans', 'Dogs' , 'God', 'Potato'};
    condition_names = {'Fluvial';'Wave';'Tidal';'All'};

    % an alternative color scheme for some plots
    c =  [0.45, 0.80, 0.69;...
          0.98, 0.40, 0.35;...
          0.55, 0.60, 0.79;...
          0.90, 0.70, 0.30]; 

    h = daboxplot(data,'xtlabels',condition_names,'colors',c,'whiskers',0,'scatter',1,'scattersize',25,'scatteralpha',0.5,'boxspacing',0.2); 
   % ylabel('Performance');
     % set x axis scale
    %if i == 2
    %    set(gca,'YScale','log');
    %end
    ylabel(label{i})
    set(gca,'FontSize',9.5);
    if i == 1
        ylim([0.6 1]);
    elseif i == 2
        ylim([0 50]);
    end
    %xl = xlim; 
    %xlim([xl(1), xl(2)+0.2]);    % make more space for the legend

end

ex_pth = "C:\Users\labuser10\OneDrive - sfu.ca\Working\SFU\Research\Modeling\conductivity_paper\dgsa_perm\connectivity_parameter_ranges"; %output destination

input_unit_name = {''; ''; '';};
find_distribution(input_var,input_idx,input_var_name,input_unit_name,ex_pth,1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make box plots for perm, slope, q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

input_idx = []; input_var = [];
input_idx = [idx_fluvial, idx_wave, idx_tidal];
input_var = [q, PERM, SLOPE];
label = {'Specific Discharge (m^3/s)'; 'Permeability (Darcy)'; 'Hydrualic gradient (m/m)'};


for i = 1:length(input_var(1,:))   
    %boxplot([input_var(idx_fluvial == 1,i), input_var(idx_wave == 1,i), input_var(idx_tidal == 1,i), input_var],)

    figure
    % data in a numreic array (+ grouping indices)
    data = NaN(length(input_var(:,i)),4);
    data(1:length(input_var(idx_fluvial == 1,i)),1) = input_var(idx_fluvial == 1,i);
    data(1:length(input_var(idx_wave == 1,i)),2) = input_var(idx_wave == 1,i);
    data(1:length(input_var(idx_tidal == 1,i)),3) = input_var(idx_tidal == 1,i);
    data(:,4) = input_var(:,i);
    %group_inx = [ones(1,30), 2.*ones(1,30), 3.*ones(1,30), 4.*ones(1,30)];

    %group_names = {'Humans', 'Dogs' , 'God', 'Potato'};
    condition_names = {'Fluvial';'Wave';'Tidal';'All'};

    % an alternative color scheme for some plots
    c =  [0.45, 0.80, 0.69;...
          0.98, 0.40, 0.35;...
          0.55, 0.60, 0.79;...
          0.90, 0.70, 0.30]; 

    h = daboxplot(data,'xtlabels',condition_names,'colors',c,'whiskers',0,'scatter',1,'scattersize',25,'scatteralpha',0.5,'boxspacing',0.2); 
   % ylabel('Performance');
     % set x axis scale
    if i == 1 || i == 2 ||  i == 3
        set(gca,'YScale','log');
    end
    ylabel(label{i})
    set(gca,'FontSize',9.5);
    xl = xlim; 
    if i == 1 
        ylim([1E-10 1E-7]);
    elseif i == 3 
        ylim([5E-5 5E-3]);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DGSA for connectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%close all

% Distance-Based Generalized Sensitivity Analysis
% Author: Celine Scheidt
% Date: August 2012
% Updated: January 2014


inputs = [discharge, wave, tide, sed_con, bath_grad, dmm, cn, shape, shore_rug];
input_names =  {'Qav';'W';'T';'Cs';'Dgrd';'Dmm';'Cn';'Shape';'SR'};
outputs = [PERM, SLOPE, C];

for ii = [5]
    
    % select model responce data
    if ii == 1
        idx = idx_fluvial;
        length(idx_fluvial(idx_fluvial == 1))
    elseif ii == 2
        idx = idx_wave;
        length(idx_wave(idx_wave == 1))
    elseif ii == 3
        idx = idx_tidal;
        length(idx_tidal(idx_tidal == 1))
    elseif ii == 4
        idx = idx_mixed;
        length(idx_mixed(idx_mixed == 1))
    else
        idx = idx_all; % all models
    end
    
    ParamsValues_before = []; ParamsNames = []; ParametersNames = [];
    ParamsValues_before = inputs(idx == 1,:);
    
    ParamsNames = input_names;
    ParametersNames = ParamsNames;

    

    jj = 1;
    for jj = 1:length(outputs(1,:))
       
        Responces = []; ParamsValues = []; 
        Responces = outputs(idx == 1,jj);
        [Responces, idx_output] = rmmissing(Responces); % remove the NaNs

        ParamsValues =  ParamsValues_before(idx_output == 0,:);
        ParametersValues = ParamsValues; % for trouble shooting
    

        %% Clustering the responce results

        %evaulate distance matrix between model responces
        d = pdist(Responces);

        %k-mediod clustering analysis
        nbclusters = 2;
        Clustering = kmedoids(d,nbclusters,20) %20 is iteration number to find optimal clustering result


        %% Evaluate the Main Effect

        %plot CDF of each parameter for each class
        %cdf_MainFactor(ParamsValues,Clustering,ParamsNames);

        %apply DGSA to evaluate the main effect
        %inputs = struct('PlotType','L1norm');
        %[PvalMainFactors,L1MainFactors,~,SensitivityMainFactorperClass] = dGSA_MainFactors(Clustering,ParamsValues,ParamsNames,inputs);

        plot_inputs = struct('PlotType','ASL');
        [PvalMainFactors,L1MainFactors,~,SensitivityMainFactorperClass] = dGSA_MainFactors_mod(Clustering,ParamsValues,ParamsNames,plot_inputs);

        %% Evaluate Sensitivity of Parameter Interactions
        %
        %define the number of bins for each conditional parameters
        NbBins = 3*ones(1,length(ParamsNames));
        NbBins(end) = 2; % for the covariance parameter

        %plot CDF of each parameter for each class
        % cdf_Interactions(ParamsValues(:,3),ParamsValues(:,1),Clustering,2,1,'T|Qav')

        %apply DGSA to evaulate conditional effect
        %inputs = struct('PlotType','L1norm','L1MainFactors',SensitivityMainFactorperClass);
        %[ALSInteractions,L1Interactions] = dGSA_Interactions(Clustering,ParamsValues,NbBins,ParamsNames,inputs);

        plot_inputs = struct('PlotType','ASL','ASLMainFactor',PvalMainFactors);
       % [ALSInteractions,L1Interactions] = dGSA_Interactions_mod(Clustering,ParamsValues,NbBins,ParamsNames,plot_inputs);
        %
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DGSA for discharge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%close all

% Distance-Based Generalized Sensitivity Analysis
% Author: Celine Scheidt
% Date: August 2012
% Updated: January 2014
%{
inputs = [cn, ncw shape, shore_rug];
input_names =  {'cn'; 'Ncw'; 'shape'; 'shore_rug'};
outputs = [q];

for ii = [5]
    
    % select model responce data
    if ii == 1
        idx = idx_fluvial;
        length(idx_fluvial(idx_fluvial == 1))
    elseif ii == 2
        idx = idx_wave;
        length(idx_wave(idx_wave == 1))
    elseif ii == 3
        idx = idx_tidal;
        length(idx_tidal(idx_tidal == 1))
    elseif ii == 4
        idx = idx_mixed;
        length(idx_mixed(idx_mixed == 1))
    else
        idx = idx_all; % all models
    end
    
    ParamsValues = [];    
    
    ParamsValues = inputs(idx == 1,:);
    ParamsNames = input_names;
        % for trouble shooting
        ParametersValues = ParamsValues;
        ParametersNames = ParamsNames;
    
    jj = 1;
    for jj = 1:length(outputs(1,:))
       
        Responces = [];
        Responces = outputs(idx == 1,jj);
        Responces = rmmissing(Responces); % remove the NaNs

        %% Clustering the responce results

        %evaulate distance matrix between model responces
        d = pdist(Responces);

        %k-mediod clustering analysis
        nbclusters = 2;
        Clustering = kmedoids(d,nbclusters,20) %20 is iteration number to find optimal clustering result


        %% Evaluate the Main Effect

        %plot CDF of each parameter for each class
        %cdf_MainFactor(ParamsValues,Clustering,ParamsNames);

        %apply DGSA to evaluate the main effect
        %inputs = struct('PlotType','L1norm');
        %[PvalMainFactors,L1MainFactors,~,SensitivityMainFactorperClass] = dGSA_MainFactors(Clustering,ParamsValues,ParamsNames,inputs);

        plot_inputs = struct('PlotType','ASL');
        [PvalMainFactors,L1MainFactors,~,SensitivityMainFactorperClass] = dGSA_MainFactors_mod(Clustering,ParamsValues,ParamsNames,plot_inputs);

        %% Evaluate Sensitivity of Parameter Interactions
        %
        %define the number of bins for each conditional parameters
        NbBins = 3*ones(1,length(ParamsNames));
        NbBins(end) = 2; % for the covariance parameter

        %plot CDF of each parameter for each class
        % cdf_Interactions(ParamsValues(:,3),ParamsValues(:,1),Clustering,2,1,'T|Qav')

        %apply DGSA to evaulate conditional effect
        %inputs = struct('PlotType','L1norm','L1MainFactors',SensitivityMainFactorperClass);
        %[ALSInteractions,L1Interactions] = dGSA_Interactions(Clustering,ParamsValues,NbBins,ParamsNames,inputs);

        plot_inputs = struct('PlotType','ASL','ASLMainFactor',PvalMainFactors);
        [ALSInteractions,L1Interactions] = dGSA_Interactions_mod(Clustering,ParamsValues,NbBins,ParamsNames,plot_inputs);
        %
    end
end
%}
