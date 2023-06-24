function [y_pdf,x_pdf,dist_type,percentile,num_samples] = fit_var_to_dist2(input,x_pdf)
% This function 1) removes outliers (if tag is on) and 2) fits appropriate
% distribution to the variable

    var = input;
    var = rmmissing(var); % remove the NaNs
    %[var,outliers] = rmoutliers(var,'quartiles'); % remove outliers from dataset
    %num_samples(i,1) = length(var);
    %num_outliers(i,1) = sum(outliers(:) ==1);
    
     if ~exist('x_pdf','var')
         % third parameter does not exist, so default it to something
          x_pdf = linspace(0,max(var),length(var)*100);
     end
    
   % 

    % count number of samples and outliers
    num_samples = length(var);
    %num_outliers(i,1) = sum(outliers(:) ==1);
        
    % find apropriate distribution
        % use Shapiro Wilks test for normality
     if swtest(var) == 0
         dist_type = {'Normal'};
     elseif swtest(log(var)) == 0
         dist_type = {'Lognormal'};
     else
         % check if distribution is exponential with a One-sample Kolmogorov-Smirnov test
         test_cdf = [var,cdf(makedist('Exponential','mu',mean(var)),var)]; % make experimental exponential function
         if kstest(var,'CDF',test_cdf) == 0
            dist_type = {'Exponential'};
         else
             disp(join("Distrubtion does not fit a normal or log normal distribution for dataset"))
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
         if isequal(dist_type{1},'Lognormal')
            var_old = var;
            var = [];
            var = var_old(var_old~=0);
         end
         pd = fitdist(var,char(dist_type)); % normal distribution for all since I have transformed the data
         y_pdf = pdf(pd,x_pdf); 
         percentile = icdf(pd,[0.10 0.25 0.375 0.5 0.625 0.75 0.90]);
     end

end

