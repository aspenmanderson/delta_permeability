% Author: Aspen Anderson
% Purpose: This code uses that data and creates figures to determine the range of plausable parameters for input into Deltf3D models
    % parameter distributions (in parameter_distributions_FROM_LIT.mat code)
    % The parameters of interst were gathered from Syvistky and Siato 2007, 
    % Neinhuis et al. 2020, and Caldwell et al. 2019
    
    % Figures of distributions supressed
    
% Last Modified: 12/2/2020

clc; close all; clear all;

set(0,'DefaultAxesFontSize', 12)
set(0,'DefaultAxesFontName','Times') 
set(0,'DefaultAxesFontWeight','Demi')
set(0,'DefaultAxesLineWidth', 1.5)
set(0,'DefaultTextFontSize', 16)
set(0,'DefaultTextFontName','Times') 
set(0,'DefaultTextFontWeight','Demi')

addpath("C:\Program Files\MATLAB\R2019b\toolbox\stats\stats\")
addpath("C:\Program Files\MATLAB\R2019b\toolbox\alchemyst-ternplot"); 
addpath("C:\Users\labuser10\OneDrive - sfu.ca\Working\SFU\Research\Modeling\delta_modeling_paper_codes\parameter_ranges_from_lit\downloaded_functions\");

% figure colors
blue = [0, 0.4470, 0.7410];
red =  [0.6350, 0.0780, 0.1840];
yellow = [0.9290, 0.6940, 0.1250];
purple = [0.4940, 0.1840, 0.5560];
green = [0.4660, 0.6740, 0.1880];
teal = [0, 0.75, 0.75];
lightblue = [0.3010, 0.7450, 0.9330];
darkgrey = [0.25, 0.25, 0.25];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

in_pth = ("C:\Users\labuser10\OneDrive - sfu.ca\Working\SFU\Research\Modeling\delta_modeling_paper_codes\parameter_ranges_from_lit\combined_data\distances_included\"); 
ex_pth = in_pth; %output destination
data_name = ("syvistky_neinhuis_caldwell_join_data_FOR_MATLAB.xlsx");

% import data
opts = detectImportOptions(strcat(in_pth,data_name));
opts.Sheet = 'data';

% select delta name
opts.SelectedVariableNames = [1]; 
%opts.DataRange = '2:11';
[delta_name] = readvars(strcat(in_pth,data_name),opts);
my_delta_ID = 1:1:length(delta_name); %deltas organized by alphabetical order
N = length(my_delta_ID);

% select input parameters
opts.SelectedVariableNames = [2 3 5 6 13 14 16 17 18 19 20 21 22 24 25]; 
[join_distance, join_point, Ad, Qav, QTide, QWave, Cs, Dmm, Wa, Wave_height, Ti, TidalAmp, Tidal_range, dgrd, Bath_slope] = readvars(strcat(in_pth,data_name),opts);

% select output parameters
opts.SelectedVariableNames = [26, 27, 29, 34, 35]; 
[Dsh, CN_Prest, ChannelSlope, Tcw, NTcw] = readvars(strcat(in_pth,data_name),opts);

% change data to NaN based on join_point (1 = within basin, 0 = point not
% within basin)
for i = 1:N
    if join_point(i) == 0
        Wave_height(i) = NaN;
        Tidal_Range(i) = NaN;
        Bath_slope(i) = NaN;
    end
end

% Average the different datasets
for i = 1:N
    if join_point(i) == 1
        Wave(i,1) = (Wa(i)+Wave_height(i))/2;
        Tide(i,1) = (Ti(i)+TidalAmp(i)+Tidal_range(i))/3;
        %Dgrd(i,1) = geomean([dgrd(i) Bath_slope(i)]);
    else
        Wave(i,1) = Wa(i);
        Tide(i,1) = (Ti(i)+TidalAmp(i))/2;
        %Dgrd(i,1) = dgrd(i);
    end
end
Dgrd = Bath_slope;

input_var = [Qav(:,1), Wave(:,1), Tide(:,1), Ad(:,1), Dgrd(:,1), Cs(:,1), str2double(Dmm(:,1))];
input_var_name = ["Qav", "Wave", "Tide", "Ad", "Bath_slope", "Cs", "Dmm"];
input_unit_name = ["m3/s", "m", "m", "km2", "m/m", "kg/m3",  "mm"];

A = invprctile(Qav,12796.4)
B = invprctile(Wa, 2.2233)
C = invprctile(Ti, 2.2749)
Rr1 = A./(A+B+C)
Rw1 = B./(A+B+C)
Rt1 = C./(A+B+C)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% View INPUT parameter distributions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

for i = 1:length(input_var(1,:))
    
    % calculate distribution
    var = [];
    var = input_var(:,i);
    var = rmmissing(var); % remove the NaNs
    %[var,outliers] = rmoutliers(var,'quartiles'); % remove outliers from dataset
    %num_samples(i,1) = length(var);
    %num_outliers(i,1) = sum(outliers(:) ==1);
    x_pdf = linspace(0,max(var),length(var)*100);

    % count number of samples and outliers
    num_samples(i,1) = length(var);
    %num_outliers(i,1) = sum(outliers(:) ==1);
        
    % find apropriate distribution
        % use Shapiro Wilks test for normality
     if swtest(var) == 0
         dist_type = {'Normal'};
         %var_transformed = var;
     elseif swtest(log(var)) == 0
         dist_type = {'Lognormal'};
         %var_transformed = log(var);
         %var_transformed = normalize(var_transformed,'range',[min(var) max(var)]);
     else
         % check if distribution is exponential with a One-sample Kolmogorov-Smirnov test
         test_cdf = [var,cdf(makedist('Exponential','mu',mean(var)),var)]; % make experimental exponential function
         if kstest(var,'CDF',test_cdf) == 0
            dist_type = {'Exponential'};
            %var_transformed = 1-exp(-0.5*var); %https://stats.stackexchange.com/questions/421559/how-to-change-exponential-distribution-into-normal-distribution
            %var_transformed = normalize(var_transformed,'range',[min(var) max(var)]);
         else
             disp(join("Distrubtion does not fit a normal or log normal distribution for dataset", i))
             dist_type = {'Kernel'};
             %var_transformed = var;
         end
     end
     

     % fit distribution
     if isequal(dist_type{1},'Kernel') 
         pd = fitdist(var,char(dist_type),'BandWidth',0.5*std(var)); % normal distribution for all since I have transformed the data
         y_pdf = pdf(pd,x_pdf);
         percentile = quantile(var,[0.10 0.25 0.375 0.50 0.625 0.75 0.90]);
     else
         pd = fitdist(var,char(dist_type)); % normal distribution for all since I have transformed the data
         y_pdf = pdf(pd,x_pdf); 
         percentile = icdf(pd,[0.10 0.25 0.375 0.5 0.625 0.75 0.90]);
     end
    
    %if i <= 3 % plot morphodynamic force parameters
        subplot(3,3,i)

        % plot dsitribution
        title(sprintf('%s', string(input_var_name(i))),'Interpreter','none')
        xlabel(sprintf('%s', string(input_unit_name(i))))

        yyaxis left
        h = histogram(var,20,'Facecolor',darkgrey); hold on
        %h = histogram(var_transformed,20, 'Facecolor', blue); hold on
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        ylabel('Frequency')
        ylim([0 max(h.Values)+1])
        ax = gca;
        ax.YAxis(1).Color = blue;

        yyaxis right
        line(x_pdf,y_pdf,'LineWidth', 2, 'Color', darkgrey); hold on %plot normal distribution
        xline(mean(var), 'LineWidth', 2, 'Color', red); %plot mean
        xline(median(var), 'LineWidth', 2, 'Color', green); %plot median
        %xline(geomean(var_transformed), 'LineWidth', 2, 'Color', yellow); %plot geomean
        ylabel(sprintf(append(char(dist_type), ' Distribution', newline,'(1/%s)'), string(input_unit_name(i))))

        %ylim([0 max([max(y_pdf), max(y_pdf_f), max(y_pdf_t)])+max([max(y_pdf), max(y_pdf_f), max(y_pdf_t)])/10])
        set(gca,'YColor','k');
        if i == 1
            legend(append('Experimental', newline, 'Distribution'),'Mean','Median','Geomean'); legend boxoff;  
        end
    
    % save parameters for table
    p_mean(i) = mean(var);
    p_med(i) = median(var);
    %p_geom(i) = geomean(var);
    p_min(i) = min(var);
    p_max(i) = max(var);
    p_iqr(i) = iqr(var); %IQR does appear to work for lognormal distributions
    p_quantile(i,:) = quantile(var,[0.10 0.25 0.375 0.50 0.625 0.75 0.90]);
    p_percentile(i,:) = percentile; 
    p_std(i) = std(var);
    %p_dist_mean(i) = pd.mu;
    %p_dist_std(i) = pd.sigma;
    p_dist(i) = dist_type;
    
end

% make table with data
T = table(input_var_name', input_unit_name', p_mean',p_med',p_min',p_max',p_std',p_dist',p_percentile(:,1),p_percentile(:,2),p_percentile(:,3),p_percentile(:,4),p_percentile(:,5),p_percentile(:,6),p_percentile(:,7), num_samples);
T.Properties.VariableNames = {'Variable','Unit','Mean','Median','Minimum','Maximum','Standard Deviation','Expirimental Distribution Type','10','25','37.5','50','62.5','75','90','Number of Samples'};
writetable(T,strcat(ex_pth,"exported_input_parameter_ranges_try7_models_new.txt"),'Delimiter',',');
writetable(T,strcat(ex_pth,"exported_input_parameter_ranges_try7_models_new.xlsx"));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Use Syviski Power Ratio to determine force distributions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fluvial, tidal, and wave are the first three parameters in p_med, p_min,
% and p_max

f = []; w = []; t = [];
cf = []; cw = []; ct = [];
Qav = []; Wa = []; Ti = [];

% select power parameteres (Syvistki)
f = [p_percentile(1,4), p_percentile(1,6), p_percentile(1,7)];
w = [0, p_percentile(2,2), p_percentile(2,4), p_percentile(2,6), p_percentile(2,7)];
t = [0, p_percentile(3,2), p_percentile(3,4), p_percentile(3,6), p_percentile(3,7)];

[cf, cw, ct] = ndgrid(f, w, t); 
combs_forces = [cf(:), cw(:), ct(:)];
Qav = combs_forces(:,1);
Wa = combs_forces(:,2);
Ti = combs_forces(:,3);

%{
n = length(Qav);
Dgrd_const = 5E-4;

nQr = []; nQw = []; nQt = [];
for i = 1:n
    nQr(i,1) = Qav(i)*Dgrd_const/((Qav(i)*Dgrd_const)+Wa(i)^2+Ti(i)^2); % normalized fluvial ratio for the Syvistki dataset
    nQw(i,1) = Wa(i)^2/((Qav(i)*Dgrd_const)+Wa(i)^2+Ti(i)^2); % normalized tidal ratio for the Syvistki dataset
    nQt(i,1) = Ti(i)^2/((Qav(i)*Dgrd_const)+Wa(i)^2+Ti(i)^2); % normalized wave ratio for the Syvistki dataset
    tot_en(i,1) = sqrt((Qav(i)*Dgrd_const)^2+Wa(i)^2+Ti(i)^2) %total energy in each system
end

% normalized total energy
for i = 1:n
    ntot_en(i,1) = tot_en(i,1)/max(tot_en)
end
%}

n = length(Qav);
for i = 1:n
    if floor(Qav(i)) == 1845
        Nr(i) = 0.5;
    elseif floor(Qav(i)) == 6768
        Nr(i) = 0.75;
    else
        Nr(i) = 0.9;
    end
    
    if Wa(i) == 0
        Nw(i) = 0;
    elseif round(Wa(i),1) == 1.0
        Nw(i) = 0.25;
    elseif round(Wa(i),1) == 1.4
        Nw(i) = 0.50;
    elseif round(Wa(i),1) == 1.9
        Nw(i) = 0.75;
    else
        Nw(i) = 0.9;
    end
    
    if Ti(i) == 0
        Nt(i) = 0;
    elseif round(Ti(i),1) == 0.4
        Nt(i) = 0.25;
    elseif round(Ti(i),1) == 1.1
        Nt(i) = 0.50;
    elseif round(Ti(i),1) == 2.1
        Nt(i) = 0.75;
    else
        Nt(i) = 0.9;
    end
    
end

Nw_Nr = Nw./Nr;

Rr = Nr./(Nr+Nw+Nt);
Rw = Nw./(Nr+Nw+Nt);
Rt = Nt./(Nr+Nw+Nt);
R_sum = Rr+Rw+Rt;




%%
figure
ternplot_mod(Rw, Rt, Rr,'sortpoints', false,'fractions', false,'.', 'markersize', 20); hold on

x = []; y = []; z = [];
[x, y] = terncoords(Rw, Rt);
z = Rr;
names = string(reshape(linspace(1,n,n),n,1));
names = strcat('F',names);
for i = 1:n 
    text(x(i)-0.01, y(i)-0.03, z(i)+0.01, names{(i)}, 'fontsize', 11);
end

% export comination in a spread sheet
T = table(Qav, Wa, Ti, Rr', Rw', Rt');
T.Properties.VariableNames = {'Fluvial discharge','Wave Height','Tidal Range','Fluvial Ratio','Wave Ratio','Tidal Ratio'};
writetable(T,strcat(ex_pth,"force_modeling_combinations_TRY7_new.txt"),'Delimiter',',');
writetable(T,strcat(ex_pth,"force_modeling_combinations_TRY7_new.xlsx"));

