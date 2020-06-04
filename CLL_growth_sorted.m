% Load in CLL untreated sorted growth data from March 20th, 2020
close all; clear all; clc;
[growth_data, T]= xlsread('../data/CLL_growth_3_20.xlsx');

%% Adjust the data
tdata = growth_data(:,1); % this is in hours
% cell numbers should be in millions
N_cells = growth_data(:,2:end).*1e6;
nsubpops = 3;
%% Make a data structure that holds each of samples
ntps = length(tdata);
CLLdata = struct('time', zeros(ntps,1), 'rawN', zeros(ntps,3), 'sample', strings(1,1),...
        'Nmean', zeros(ntps,1), 'Nstd', zeros(ntps,1));
    
% Go through and add each sample
sampsnames = {'TP0', 'TP0 CD18+', 'TP0 CXCR4+'};
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
%% Fit to a single logistig growth model to set the carrying capacity going forward
Nall = [];
for i = 1:nsubpops
    Nall = horzcat(Nall,CLLdata(i).rawN);
end
Nallmean = mean(Nall,2);
Nallstd = std(Nall,0,2);
Nallstd(1) = 1e4;
sigma = 1e4;
[punt, singexpmodel] = fit_untreated(Nallmean,tdata, sigma);
gglobal = punt(1)
Kglob = punt(2)
figure;
errorbar(tdata, Nallmean, 1.96*Nallstd, 'b*-')
hold on
plot(tdata, singexpmodel, 'k-')
xlabel('time(hours)')
ylabel('Global fit for best K')
set(gca,'FontSize',20,'LineWidth',1.5)
%% Set K and fit for g
%K = 15e7;
for i = 1:nsubpops
    sigmavec = CLLdata(i).Nstd;
    sigma = 1e4; % add a little uncertainty to initial condition
    [p, modelfitg] = fit_logistic(CLLdata(i).Nmean, tdata, sigma, Kglob);
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
        xlim([tdata(1) tdata(end)])
        legend(sampsnames(i), 'Location', 'NorthWest')
        legend boxoff
        ylim([0 1e8])
        
end

%% Fit individually for k and g

for i = 1:nsubpops
    sigmavec = CLLdata(i).Nstd;
    sigma = 1e4; % add a little uncertainty to initial condition
    [p, modelfit] = fit_untreated(CLLdata(i).Nmean, tdata, sigma);
    CLLdata(i).gandK = p;
    CLLdata(i).modelfitwithk = modelfit;
    CLLdata(i).doublingtime = log(2)./p(1)
end

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
        title(['g=', num2str(round(CLLdata(i).gandK(1),4)), ', t_{double}=', num2str(CLLdata(i).doublingtime),' hrs'])%,'K=', num2str(CLLdata(i).gandK(2))])
        set(gca,'FontSize',20,'LineWidth',1.5)
        xlim([tdata(1) tdata(end)])
        legend(sampsnames(i), 'Location', 'NorthWest')
        legend boxoff
        ylim([0 1e8])
        
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
% CXCR4 vs CD18 stat significant?
[h,p1] = ttest2(CLLdata(2).gset, CLLdata(3).gset)
% Are TP0 bulk and CD18 stat significant different?
[h,p2] = ttest2(CLLdata(1).gset, CLLdata(2).gset)
% Are TP0 bulk and CXCR4 stat significant differetn?
[h,p3] = ttest2(CLLdata(1).gset, CLLdata(3).gset)

%% figure to emulate Cathy's
figure;
errorbar(CLLdata(1).time, CLLdata(1).Nmean,1.96*CLLdata(1).Nstd, '*', 'color', CLLdata(1).color, 'LineWidth', 2)
 hold on
errorbar(CLLdata(2).time, CLLdata(2).Nmean,1.96*CLLdata(2).Nstd, '*', 'color', CLLdata(2).color, 'LineWidth', 2)
errorbar(CLLdata(3).time, CLLdata(3).Nmean,1.96*CLLdata(3).Nstd, '*', 'color', CLLdata(3).color, 'LineWidth', 2)
plot(CLLdata(1).time, CLLdata(1).modelfitwithk, '-', 'color', CLLdata(1).color, 'LineWidth', 2) 
plot(CLLdata(2).time, CLLdata(2).modelfitwithk, '-', 'color', CLLdata(2).color, 'LineWidth', 2) 
plot(CLLdata(3).time, CLLdata(3).modelfitwithk, '-', 'color', CLLdata(3).color, 'LineWidth', 2) 
xlabel('time (hours)')
ylabel('N(t)')
title( 'TP0')
%title(['g=', num2str(round(CLLdata(i).gandK(1),2)), ', t_{double}=', num2str(round(CLLdata(i).doublingtime),3),' hrs'])%,'K=', num2str(CLLdata(i).gandK(2))])
set(gca,'FontSize',20,'LineWidth',1.5)
xlim([tdata(1) tdata(end)])
legend('TP0 t_{double} = 24 hrs', 'TP0 CD18^{hi} t_{double} = 31 hrs', 'TP0 CXCR4^{hi} t_{double} = 19 hrs', 'Location', 'NorthWest')
legend boxoff
ylim([0 1.2e8])