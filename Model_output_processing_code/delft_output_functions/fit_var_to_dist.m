function [y_pdf,x_pdf,dist_type] = fit_var_to_dist(var)
% This function 1) removes outliers (if tag is on) and 2) fits appropriate
% distribution to the variable

    var = rmmissing(var); % remove the NaNs
    var = nonzeros(var);
    [var,outliers] = rmoutliers(var,'quartiles'); % remove outliers from dataset
    %num_samples(i,1) = length(var);
    %num_outliers(i,1) = sum(outliers(:) ==1);
    x_pdf = linspace(0,max(var),length(var)*100);


    % find apropriate distribution for literature variable
        % use Shapiro Wilks test for normality
     if swtest(var) == 0
         dist_type = {'Normal'};
         %[var,outliers] = rmoutliers(var,'percentiles',[5 95]);
     elseif swtest(log(var)) == 0
         dist_type = {'Lognormal'};
         %[var,outliers] = rmoutliers(log(var),'percentiles',[5 95]);
         %var = exp(var);
     else
         % check if distribution is exponential with a One-sample Kolmogorov-Smirnov test
         test_cdf = [var,cdf(makedist('Exponential','mu',mean(var)),var)]; % make experimental exponential function
         if kstest(var,'CDF',test_cdf) == 0
            dist_type = {'Exponential'};
            %[var,outliers] = rmoutliers(var,'percentiles',[5 95]);
         else
             dist_type = {'Kernel'};
             [var,outliers] = rmoutliers(var,'percentiles',[5 95]);
         end
     end

     % count number of samples and outliers
     num_samples = length(var);
     num_outliers = sum(outliers(:) ==1);

     %pd = fitdist(var,char(dist_type)); % normal distribution for all since I have transformed the data
     %y_pdf = pdf(pd,x_pdf);
     if isequal(dist_type{1},'Kernel') 
         pd = fitdist(var,char(dist_type),'BandWidth',0.5*std(var)); % normal distribution for all since I have transformed the data
         y_pdf = pdf(pd,x_pdf);
     else
         pd = fitdist(var,char(dist_type)); % normal distribution for all since I have transformed the data
         y_pdf = pdf(pd,x_pdf); 
     end

end

