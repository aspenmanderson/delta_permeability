function [] = find_distribution(input_var,idx,input_var_name,input_unit_name,ex_pth,type)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    % figure colors
    blue = [0, 0.4470, 0.7410];
    red =  [0.6350, 0.0780, 0.1840];
    yellow = [0.9290, 0.6940, 0.1250];
    purple = [0.4940, 0.1840, 0.5560];
    green = [0.4660, 0.6740, 0.1880];
    teal = [0, 0.75, 0.75];
    lightblue = [0.3010, 0.7450, 0.9330];
    darkgrey = [0.25, 0.25, 0.25];
%%
    %Perm, AD_hp*100, slope, CN, shape, shorerug
    %log axis = perm & slope
    %{
    figure
    %label = {'Fluvial Discharge (m^3/s)'; 'Wave Height (m)'; 'Tidal Range (m)'; 'Bathymetric Gradient (m/km)'; 'Sediment Concentration (kg/m^3)'; 'Median Grain Size (mm)'};
    label = input_var_name;

    for i = 1:length(input_var(1,:))
        
        binNum = 40;
        threshold = [0 99];
        
        % make bins
        if i == 1 || i == 2 || i == 3
            binRange = logspace(log10(min(input_var(idx(:,5) == 1,i))),log10(max(input_var(idx(:,5) == 1,i))),binNum);
            binRange = sort(binRange);
            hcy = histcounts(input_var(idx(:,5) == 1,i),[binRange Inf]);
        else 
            binRange = 0:abs(max(input_var(idx(:,5) == 1,i))/binNum):max(input_var(idx(:,5) == 1,i));
            binRange = sort(binRange);
            hcy = histcounts(input_var(idx(:,5) == 1,i),[binRange Inf]);
        end
        
        x_pdf = []; y_pdf = []; 
        y_pdf_f = []; y_pdf_w = []; y_pdf_t = []; y_pdf_m = [];
        % create distributions
        if i == 1 || i == 2 || i == 3   
            x_pdf = logspace(log10(binRange(1)),log10(binRange(end)),length(input_var(idx(:,5) == 1,i))*1000);
        else          
            x_pdf = linspace(binRange(1),binRange(end),length(input_var(idx(:,5) == 1,i))*100);     
        end     
        [y_pdf_f] = fit_var_to_dist2(input_var(idx(:,1) == 1,i),x_pdf);  % fluvial          
        [y_pdf_w] = fit_var_to_dist2(input_var(idx(:,2) == 1,i),x_pdf);  % wave           
        [y_pdf_t] = fit_var_to_dist2(input_var(idx(:,3) == 1,i),x_pdf); % tidal 
        [y_pdf_m] = fit_var_to_dist2(input_var(idx(:,4) == 1,i),x_pdf); % mixed
        [y_pdf,x_pdf,dist_type,percentile] = fit_var_to_dist2(input_var(idx(:,5) == 1,i),x_pdf); % all
    

        % make figure
        subplot(3,3,i)
        hold on;
        xlabel(label(i))
        ylabel('Proportion of Data')

        % plot histograms
        histogram(input_var(idx(:,5) == 1,i),binRange,'Normalization','probability','FaceAlpha',0.7,'FaceColor',darkgrey); hold on % create the plot 

        % plot input parameter dots
       % scatter(inv,zeros(length(inv),1),40,'y','filled', 'MarkerEdgeColor','k');

        % y axis
        %ylim([0 max(hcy/length(input_var))+max(hcy/length(input_var))*0.1])

        % x axis
        if i == 1 || i == 2 || i == 3 
            set(gca,'XScale','log');
        end
        edge_low = binRange(2)-binRange(1);
        edge_high = binRange(end)-binRange(end-1);
        xlim([min(input_var(idx(:,5) == 1,i))-edge_low max(input_var(idx(:,5) == 1,i))+edge_high]);



        % plot distributions
        yyaxis right
        %
        if i == 2 || i == 3 
            line(x_pdf,y_pdf_f,'LineWidth', 2, 'Color', red); hold on %plot normal distribution
            line(x_pdf,y_pdf_w,'LineWidth', 2, 'Color', blue); 
            line(x_pdf,y_pdf_t,'LineWidth', 2, 'Color', yellow); 
            line(x_pdf,y_pdf_m,'LineWidth', 2, 'Color', green); 
        else
            line(x_pdf,y_pdf_f,'LineWidth', 2, 'Color', red); hold on %plot normal distribution
            line(x_pdf,y_pdf_w,'LineWidth', 2, 'Color', blue); 
            line(x_pdf,y_pdf_t,'LineWidth', 2, 'Color', yellow); 
            line(x_pdf,y_pdf_m,'LineWidth', 2, 'Color', green); 
        end
        %
      
        ylabel(sprintf('Fit Distribution: (1/%s)', string(input_unit_name(i)))) 

        if i == 2 || i == 3  
            set(gca,'XScale','log');
        end

        ax = gca;
        ax.YAxis(1).Color = 'k';
        ax.YAxis(2).Color = 'k';
        %

       if i == 2
            legend('All models','Fluvial','Wave','Tidal','Mixed'); legend boxoff;
       end
    
%{
    % save parameters for table
    p_mean(i) = mean(var);
    p_med(i) = median(var);
    %p_geom(i) = geomean(var);
    p_min(i) = min(var);
    p_max(i) = max(var);
    p_percentile(i,:) = percentile; 
    p_std(i) = std(var);
    %p_dist_mean(i) = pd.mu;
    %p_dist_std(i) = pd.sigma;
    p_dist(i) = dist_type;
    
    % make table with data
    T = table(input_var_name', input_unit_name', p_mean',p_med',p_min',p_max',p_std',p_dist',p_percentile(:,1),p_percentile(:,2),p_percentile(:,3),p_percentile(:,4),p_percentile(:,5),p_percentile(:,6),p_percentile(:,7), num_samples);
    T.Properties.VariableNames = {'Variable','Unit','Mean','Median','Minimum','Maximum','Standard Deviation','Expirimental Distribution Type','10','25','35.75','50','62.25','75','90','Number of Samples'};
    writetable(T,strcat(ex_pth,".txt"),'Delimiter',',');
    writetable(T,strcat(ex_pth,".xlsx"));
    %}
    % spread sheet for all delta models at last time step
    end
    %}
    
    if type == 1
        table_name = strcat(ex_pth,".xlsx");
        for i = 1:length(input_var(1,:))
            for j = 1:length(idx(1,:))
                var = input_var(idx(:,j) == 1,i);
                data_table = [nanmean(var), median(var,'omitnan'), geomean(abs(var),'omitnan'), min(var) max(var) std(var,'omitnan') prctile(var,10) prctile(var,25) prctile(var,50) prctile(var,75) prctile(var,90)];% dist_type percentile(:,1) percentile(:,2) percentile(:,3) percentile(:,4) percentile(:,5) percentile(:,6) percentile(:,7)]; 
                xlswrite(table_name,data_table,string(input_var_name(i)),num2str(j)); % all
            end
        end
    end
    if type == 2
        table_name = strcat(ex_pth,".xlsx");
        for i = 1:length(idx(1,:))
            for j = 1:length(input_var(1,:))
                var = input_var(idx(:,i) == 1,:);
                data_table = [nanmean(var), median(var,'omitnan'), geomean(abs(var),'omitnan'), min(var) max(var) std(var,'omitnan') prctile(var,10) prctile(var,25) prctile(var,50) prctile(var,75) prctile(var,90)];% dist_type percentile(:,1) percentile(:,2) percentile(:,3) percentile(:,4) percentile(:,5) percentile(:,6) percentile(:,7)]; 
                xlswrite(table_name,data_table,string(input_var_name(i)),num2str(j)); % all
            end
        end
    end


   
end

