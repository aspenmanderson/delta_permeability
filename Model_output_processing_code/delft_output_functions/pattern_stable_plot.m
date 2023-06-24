function [] = pattern_stable_plot(var, type, morpho_time)
%UNTITLED3 Summary of this function goes here
%  Not finished
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Has the delta reached pattern stable?
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % indicators of pattern stable 
        % - steady incres in AD, shape?,
        % number of islands, NTcw, CN, Dgrd, shoreline rugosity
    
    %{     
    % calucating how much the moving mean is changing by
    for i = 1:length(Tcw)
        if Tcw(i) == 0
            NTcw(i) = 0;
        else
            NTcw(i) = Tcw(i)/CN(i);
        end
    end
    
    
    figure
    subplot(3,3,1); 
    pattern_stable_plot(AD,'linear',morpho_time)
    legend('AD'); legend boxoff

    subplot(3,3,2)
    pattern_stable_plot(L,'linear',morpho_time)
    legend('L'); legend boxoff

    subplot(3,3,3)
    pattern_stable_plot(shape,'linear',morpho_time)
    legend('shape'); legend boxoff

    subplot(3,3,4)
    pattern_stable_plot(CN,'mean',morpho_time)
    legend('CN'); legend boxoff

    subplot(3,3,5)
    pattern_stable_plot(num_islands,'mean',morpho_time)
    legend('num_islands'); legend boxoff

    subplot(3,3,6)
    pattern_stable_plot(NTcw,'mean',morpho_time)
    legend('NTcw'); legend boxoff

    subplot(3,3,7)
    pattern_stable_plot(shore_rug,'mean',morpho_time)
    legend('shore_rug'); legend boxoff;

    subplot(3,3,8)
    pattern_stable_plot(Dgrd,'mean',morpho_time)
    legend('Dgrd'); legend boxoff;

    subplot(3,3,9)
    pattern_stable_plot(Dsh,'linear',morpho_time)
    legend('Dsh'); legend boxoff;
    %}
        
    red =  [0.80, 0.00, 0.20];
    green = [0.4660, 0.6740, 0.1880];
    purple = [0.4940, 0.1840, 0.5560];
    
    %var = diff(var)./var(1:end-1,:);
    
    if strcmp(type,'linear')
        scatter(morpho_time/morpho_time(end),var,'.'); hold on 
        ipt = findchangepts(var,'Statistic','linear','MinThreshold',0.25*(var(end)-var(1)));
        xline(morpho_time(ipt(end))/morpho_time(end), 'LineWidth', 1, 'Color', red);

        y = nan(size(var));
        K = length(ipt);
        nseg = K+1;
        istart = [1; ipt(:)];
        istop = [ipt(:)-1; length(var)];
        for s = 1:nseg
                ix = (istart(s):istop(s))';
                y(ix) = polyval(polyfit(ix,var(ix)',1),ix);
                plot(morpho_time(ix)/morpho_time(end),y(ix),'Color',purple)
        end
        
    elseif strcmp(type,'mean')
        scatter(morpho_time/morpho_time(end),var,'.'); hold on 
        ipt = findchangepts(var,'Statistic','mean');
        xline(morpho_time(ipt(end))/morpho_time(end), 'LineWidth', 1, 'Color', red);

        y = nan(size(var));
        K = length(ipt);
        nseg = K+1;
        istart = [1; ipt(:)];
        istop = [ipt(:)-1; length(var)];
        for s = 1:nseg
                ix = (istart(s):istop(s))';
                y(ix) = mean(var(ix));
                plot(morpho_time(ix)/morpho_time(end),y(ix),'Color',purple)
        end
        
    end

end

