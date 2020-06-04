% Load in CD18 and CXCR4 sorted growth data from March 20th, 2020
close all; clear all; clc;
[CD18_data, T1]= xlsread('../data/CD18_growth_data.xlsx');
[CXCR4_data, T2] = xlsread('../data/CXCR4_growth_data.xlsx');
%% Adjust the data
tCD18 = CD18_data(:,1); % this is in hours
tCXCR4 = CXCR4_data(:,1); 
nsubpops = 4;
%% Make a data structure that holds each of samples
ntps = length(tCD18);
CLLdata = struct('time', zeros(ntps,1), 'rawN', zeros(ntps,3), 'sample', strings(1,1),...
        'Nmean', zeros(ntps,1), 'Nstd', zeros(ntps,1));
    
% Go through and add each sample
sampsnames = {'TP0 CD18+', 'TP1 CD18+', 'TP0 CXCR4+', 'TP1 CXCR4+'};
colorsets = varycolor(nsubpops); % give each sample a unique color (for plotting)
CLLdata(1).time = tCD18;
CLLdata(1).rawN = CD18_data(:,2:4).*1e6;
CLLdata(2).time = tCD18;
CLLdata(2).rawN = CD18_data(:,5:7).*1e6;
CLLdata(3).time = tCXCR4;
CLLdata(3).rawN = CXCR4_data(:,2:4).*1e6;
CLLdata(4).time = tCXCR4;
CLLdata(4).rawN = CXCR4_data(:,5:7).*1e6;
for i = 1:nsubpops
    CLLdata(i).sample = sampsnames(i);
    CLLdata(i).Nmean = nanmean(CLLdata(i).rawN,2);
    CLLdata(i).Nstd = nanstd(CLLdata(i).rawN, 0, 2);
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
        xlim([tCD18(1) tCD18(end)])
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
        xlim([tCD18(1) tCD18(end)])
        legend(sampsnames, 'Location', 'NorthWest')
        legend boxoff
end
%%
figure;
for i = 1:nsubpops
    for j = 1:3
    plot(CLLdata(i).time, CLLdata(i).rawN(:,j), '*-', 'color', CLLdata(i).color, 'LineWidth', 2) 
    hold on
    end
         %text(CLLdata(i).time(end-2), CLLdata(i).Nmean(end-2), [(CLLdata(i).sample)])
         %plot(CLLdata(i).time, CLLdata(i).Nmean + 1.96*CLLdata(i).Nstd, 'color', CLLdata(i).color)
         %plot(CLLdata(i).time, CLLdata(i).Nmean - 1.96*CLLdata(i).Nstd, 'color', CLLdata(i).color)
        xlabel('time (hours)')
        ylabel('N(t)')
        title('N(t) all data colored by sample')
        set(gca,'FontSize',20,'LineWidth',1.5)
        xlim([tCD18(1) tCD18(end)])
        %legend(sampsnames, 'Location', 'NorthWest')
        %legend boxoff
end

%% Set K and fit for g
K = 15e7;
for i = 1:nsubpops
    sigmavec = CLLdata(i).Nstd;
    sigma = 1e4; % add a little uncertainty to initial condition
    [p, modelfitg] = fit_logistic(CLLdata(i).Nmean, CLLdata(i).time, sigma, K);
    CLLdata(i).g = p(1);
    CLLdata(i).modelfitg = modelfitg;
end
figure;
for i = 1:nsubpops
    subplot(1, nsubpops, i)
      errorbar(CLLdata(i).time, CLLdata(i).Nmean,1.96*CLLdata(i).Nstd, '*', 'color', CLLdata(i).color, 'LineWidth', 2)
         hold on
         plot(CLLdata(i).time, CLLdata(i).modelfitg, '-', 'color', CLLdata(i).color, 'LineWidth', 2) 
         %text(CLLdata(i).time(end-2), CLLdata(i).modelfit(end-2), ['g=', num2str(CLLdata(i).gfit)])
         %plot(CLLdata(i).time, CLLdata(i).Nmean + 1.96*CLLdata(i).Nstd, 'color', CLLdata(i).color)
         %plot(CLLdata(i).time, CLLdata(i).Nmean - 1.96*CLLdata(i).Nstd, 'color', CLLdata(i).color)
        xlabel('time (hours)')
        ylabel('N(t)')
        title(['g=', num2str(CLLdata(i).g)])
        set(gca,'FontSize',20,'LineWidth',1.5)
        xlim([tCD18(1) tCD18(end)])
        legend(sampsnames(i), 'Location', 'NorthWest')
        legend boxoff
        ylim([0 1e8])
        
end

%% Fit individually for k and g

for i = 1:nsubpops
    sigmavec = CLLdata(i).Nstd;
    sigma = 1e4; % add a little uncertainty to initial condition
    [p, modelfit] = fit_untreated(CLLdata(i).Nmean, CLLdata(i).time, sigma);
    CLLdata(i).gandK = p;
    CLLdata(i).modelfitwithk = modelfit;
    CLLdata(i).doublingtime = log(2)./p(1);
end

%% Fit each replicate individually for g and k
for i = 1:nsubpops
    sigma = 1e4;
    gset = [];
    for j = 1:3
        Ni = [];
        tvec = [];
        Ni = CLLdata(i).rawN(~isnan(CLLdata(i).rawN(:,j)),j);
        tvec = CLLdata(i).time(~isnan(CLLdata(i).rawN(:,j)));
    [pi, modelfiti] = fit_untreated(Ni,tvec, sigma);
    gset(j) = pi(1)
    end
    CLLdata(i).gset = gset;
end
[hCD18,pCD18] = ttest2(CLLdata(1).gset, CLLdata(2).gset)
[hCXCR4,pCXCR4] = ttest2(CLLdata(3).gset, CLLdata(4).gset)

%%
figure;
for i = 1:nsubpops
    subplot(1, nsubpops, i)
      errorbar(CLLdata(i).time, CLLdata(i).Nmean,1.96*CLLdata(i).Nstd, '*', 'color', CLLdata(i).color, 'LineWidth', 2)
         hold on
         plot(CLLdata(i).time, CLLdata(i).modelfitwithk, '-', 'color', CLLdata(i).color, 'LineWidth', 2) 
         %text(CLLdata(i).time(end-2), CLLdata(i).modelfit(end-2), ['g=', num2str(CLLdata(i).gfit)])
         %plot(CLLdata(i).time, CLLdata(i).Nmean + 1.96*CLLdata(i).Nstd, 'color', CLLdata(i).color)
         %plot(CLLdata(i).time, CLLdata(i).Nmean - 1.96*CLLdata(i).Nstd, 'color', CLLdata(i).color)
        xlabel('time (hours)')
        ylabel('N(t)')
        title(['g=', num2str(round(CLLdata(i).gandK(1),2)), ', t_{double}=', num2str(round(CLLdata(i).doublingtime),3),' hrs'])%,'K=', num2str(CLLdata(i).gandK(2))])
        set(gca,'FontSize',20,'LineWidth',1.5)
        xlim([tCD18(1) tCD18(end)])
        legend(sampsnames(i), 'Location', 'NorthWest')
        legend boxoff
        ylim([0 1e8])
        
end
%% figure to emulate Cathy's
figure;
errorbar(CLLdata(1).time, CLLdata(1).Nmean,1.96*CLLdata(1).Nstd, '*', 'color', 'b', 'LineWidth', 2)
 hold on
errorbar(CLLdata(2).time, CLLdata(2).Nmean,1.96*CLLdata(2).Nstd, '*', 'color', [0.5, 0 , 0.5], 'LineWidth', 2)
plot(CLLdata(1).time, CLLdata(1).modelfitwithk, '-', 'color', 'b', 'LineWidth', 2) 
plot(CLLdata(2).time, CLLdata(2).modelfitwithk, '-', 'color', [0.5, 0, 0.5], 'LineWidth', 2) 
xlabel('time (hours)')
ylabel('N(t)')
title( 'CD18^{hi}')
%title(['g=', num2str(round(CLLdata(i).gandK(1),2)), ', t_{double}=', num2str(round(CLLdata(i).doublingtime),3),' hrs'])%,'K=', num2str(CLLdata(i).gandK(2))])
set(gca,'FontSize',20,'LineWidth',1.5)
xlim([tCD18(1) tCD18(end)])
legend('TP0 CD18^{hi} t_{double} = 31 hrs', 'TP1 CD18^{hi} t_{double} = 18 hrs', 'Location', 'NorthWest')
legend boxoff
ylim([0 8e7])

figure;
errorbar(CLLdata(3).time, CLLdata(3).Nmean,1.96*CLLdata(3).Nstd, '*', 'color', 'r', 'LineWidth', 2)
 hold on
errorbar(CLLdata(4).time, CLLdata(4).Nmean,1.96*CLLdata(4).Nstd, '*', 'color', [ 0.9100 0.4100 0.1700], 'LineWidth', 2)
plot(CLLdata(3).time, CLLdata(3).modelfitwithk, '-', 'color', 'r', 'LineWidth', 2) 
plot(CLLdata(4).time, CLLdata(4).modelfitwithk, '-', 'color', [ 0.9100 0.4100 0.1700], 'LineWidth', 2) 
xlabel('time (hours)')
ylabel('N(t)')
title( 'CXCR4^{hi}')
%title(['g=', num2str(round(CLLdata(i).gandK(1),2)), ', t_{double}=', num2str(round(CLLdata(i).doublingtime),3),' hrs'])%,'K=', num2str(CLLdata(i).gandK(2))])
set(gca,'FontSize',20,'LineWidth',1.5)
xlim([tCD18(1) tCD18(end)])
legend('TP0 CXCR4^{hi} t_{double} = 18 hrs', 'TP1 CXCR4^{hi} t_{double} = 17 hrs', 'Location', 'NorthWest')
legend boxoff
ylim([0 9e7])
      