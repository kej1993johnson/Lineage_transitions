% Load in the TP0 treated data set
% p-value analysis

close all; clear all; clc;
[growth_data, T]= xlsread('../data/TP0_sorted_treated_4_9_2020.xlsx');

%% Adjust the data
tdata = growth_data(:,1); % this is in days
% cell numbers should be in millions
N_cells = growth_data(:,2:end).*1e6;
nsubpops = 3;
%% Make a data structure that holds each of samples
ntps = length(tdata);
CLLdata = struct('time', zeros(ntps,1), 'rawN', zeros(ntps,3), 'sample', strings(1,1),...
        'Nmean', zeros(ntps,1), 'Nstd', zeros(ntps,1));
    
% Go through and add each sample
sampsnames = {'TP0 FM', 'TP0 CD18+ FM', 'TP0 CXCR4+FM'};
colorsets = [0 0 0; 0 0 1; 1 0 0]; % give each sample a unique color (for plotting)
for i = 1:nsubpops
    CLLdata(i).time = tdata;
    CLLdata(i).sample = sampsnames(i);
    CLLdata(i).rawN = N_cells(:, i*3-2:i*3);
    CLLdata(i).Nmean = mean(CLLdata(i).rawN,2);
    CLLdata(i).Nstd = std(CLLdata(i).rawN, 0, 2);
    CLLdata(i).color =colorsets(i,:);
end

%% Make a plot of the mean and stdev of data
% Plot the mean of the data
figure;
for i = 1:nsubpops
    plot(CLLdata(i).time, CLLdata(i).Nmean,'*-', 'color', CLLdata(i).color, 'LineWidth', 2)
         hold on
         %text(CLLdata(i).time(end-2), CLLdata(i).Nmean(end-2), [(CLLdata(i).sample)])
         %plot(CLLdata(i).time, CLLdata(i).Nmean + 1.96*CLLdata(i).Nstd, 'color', CLLdata(i).color)
         %plot(CLLdata(i).time, CLLdata(i).Nmean - 1.96*CLLdata(i).Nstd, 'color', CLLdata(i).color)
        xlabel('time (hours)')
        ylabel('N(t)')
        title('<N(t)> for each sample')
        set(gca,'FontSize',20,'LineWidth',1.5)
        xlim([tdata(1) tdata(end)])
        legend(sampsnames, 'Location', 'NorthWest')
        legend boxoff
end
figure;
for i = 1:nsubpops
      errorbar(CLLdata(i).time, CLLdata(i).Nmean,1.96*CLLdata(i).Nstd, '*-', 'color', CLLdata(i).color, 'LineWidth', 2)
         hold on
         %text(CLLdata(i).time(end-2), CLLdata(i).Nmean(end-2), [(CLLdata(i).sample)])
         %plot(CLLdata(i).time, CLLdata(i).Nmean + 1.96*CLLdata(i).Nstd, 'color', CLLdata(i).color)
         %plot(CLLdata(i).time, CLLdata(i).Nmean - 1.96*CLLdata(i).Nstd, 'color', CLLdata(i).color)
        xlabel('time (hours)')
        ylabel('N(t)')
        title('N(t) for each sample')
        set(gca,'FontSize',20,'LineWidth',1.5)
        xlim([tdata(1) tdata(end)])
        legend(sampsnames, 'Location', 'NorthWest')
        legend boxoff
end
%%
figure;
for i = 1:nsubpops
    for j = 1:3
    plot(CLLdata(i).time, CLLdata(i).rawN(:,j), '*-', 'color', CLLdata(i).color, 'LineWidth', 2)
    end  
    hold on
         %text(CLLdata(i).time(end-2), CLLdata(i).Nmean(end-2), [(CLLdata(i).sample)])
         %plot(CLLdata(i).time, CLLdata(i).Nmean + 1.96*CLLdata(i).Nstd, 'color', CLLdata(i).color)
         %plot(CLLdata(i).time, CLLdata(i).Nmean - 1.96*CLLdata(i).Nstd, 'color', CLLdata(i).color)
        xlabel('time (hours)')
        ylabel('N(t)')
        title('N(t) all data colored by sample')
        set(gca,'FontSize',20,'LineWidth',1.5)
        xlim([tdata(1) tdata(end)])
        %legend(sampsnames, 'Location', 'NorthWest')
        %legend boxoff
end

%% Set K and fit for g
K = 5e6;
for i = 1:nsubpops
    sigmavec = CLLdata(i).Nstd;
    sigma = 1e4; % add a little uncertainty to initial condition
    [p, modelfitg] = fit_logistic(CLLdata(i).Nmean(4:end), tdata(4:end), sigma, K);
    CLLdata(i).g = p(1);
    CLLdata(i).modelfitg = modelfitg;
end
figure;
for i = 1:nsubpops
    subplot(1, nsubpops, i)
      errorbar(CLLdata(i).time(4:end), CLLdata(i).Nmean(4:end),1.96*CLLdata(i).Nstd(4:end), '*', 'color', CLLdata(i).color, 'LineWidth', 2)
         hold on
         plot(CLLdata(i).time(4:end), CLLdata(i).modelfitg, '-', 'color', CLLdata(i).color, 'LineWidth', 2) 
         %text(CLLdata(i).time(end-2), CLLdata(i).modelfit(end-2), ['g=', num2str(CLLdata(i).gfit)])
         %plot(CLLdata(i).time, CLLdata(i).Nmean + 1.96*CLLdata(i).Nstd, 'color', CLLdata(i).color)
         %plot(CLLdata(i).time, CLLdata(i).Nmean - 1.96*CLLdata(i).Nstd, 'color', CLLdata(i).color)
        xlabel('time (hours)')
        ylabel('N(t)')
        title(['g=', num2str(CLLdata(i).g)])
        set(gca,'FontSize',20,'LineWidth',1.5)
        xlim([tdata(4) tdata(end)])
        legend(sampsnames(i), 'Location', 'NorthWest')
        legend boxoff
        %ylim([0 1e8])
        
end


%% Fit each replicate individually for g and k
figure;
for i = 1:nsubpops
    sigma = 1e4;
    gset = [];
    subplot(1, 3, i)
    for j = 1:3
        Ni = [];
        tvec = [];
        Ni = CLLdata(i).rawN(~isnan(CLLdata(i).rawN(:,j)),j);
        Ni = Ni(4:end);
        tvec = CLLdata(i).time(~isnan(CLLdata(i).rawN(:,j)));
        tvec = tvec(4:end)-8;
    [pi, modelfiti] = fit_logistic(Ni,tvec, sigma, K);
    gset(j) = pi(1)
    plot(tvec+8, Ni, '*', 'color', CLLdata(i).color)
    hold on
    plot(tvec+8, modelfiti,'color', CLLdata(i).color)
    end
    CLLdata(i).gset = gset;
    xlabel('time (days)')
    ylabel('N(t)')
    title([(CLLdata(i).sample)])
    set(gca,'FontSize',20,'LineWidth',1.5)
    ylim([0 5e6])
   
end
% CXCR4 vs CD18 stat significant?
[h,p1] = ttest2(CLLdata(2).gset, CLLdata(3).gset)
% Are TP0 bulk and CD18 stat significant different?
[h,p2] = ttest2(CLLdata(1).gset, CLLdata(2).gset)
% Are TP0 bulk and CXCR4 stat significant differetn?
[h,p3] = ttest2(CLLdata(1).gset, CLLdata(3).gset)

