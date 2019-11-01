% Growth analysis of isolated subpopulations of CLL cells

% This script loads in the CLL growth rate data of the followign isolated
% populations:

% TPO (8 days post thaw, t=0 for scRNAseq)
% TP0 CD18+ (cells from this TP0 population sorted into by high CD18 
% expression levels)
% TP0 CXCR4 (cells from TP0 sorted by high CXCR4 levels (thought to be high
% tolerance)
% FM1 (cells post fluderamine treatment from sample 1)
% FM1 CD18 (cells post fluderamine treatment sorted for high CD18)
% FM1 CXCR4 (cells from fluderamine treatment sorted for high CXCR4)
% FM7 another sample of cells post fluderamine treatment ?
close all; clear all; clc;
[growth_data, T]= xlsread('../data/CLL_growth_data.xlsx');
%% Adjust the data
tdata = growth_data(:,1); % this is in hours
% cell numbers should be in millions
N_cells = growth_data(:,2:end).*1e6;
nsubpops = 7;
%% Make a data structure that holds each of samples
ntps = length(tdata);
CLLdata = struct('time', zeros(ntps,1), 'rawN', zeros(ntps,3), 'sample', strings(1,1),...
        'Nmean', zeros(ntps,1), 'Nstd', zeros(ntps,1));
    
% Go through and add each sample
sampsnames = {'TP0', 'TP0 CD18+', 'TP0 CXCR4+', 'FM1', 'FM1 CD18+', 'FM1 CXR4+', 'FM7'};
colorsets = varycolor(nsubpops); % give each sample a unique color (for plotting)
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
        title(['g=', num2str(round(CLLdata(i).gandK(1),4))])%,'K=', num2str(CLLdata(i).gandK(2))])
        set(gca,'FontSize',20,'LineWidth',1.5)
        xlim([tdata(1) tdata(end)])
        legend(sampsnames(i), 'Location', 'NorthWest')
        legend boxoff
        ylim([0 1e8])
        
end


% Save the carrying capacity for the untreated control to be loaded in
Kunt = CLLdata(1).gandK(2);
save('../out/KuntCLLT25.mat', 'Kunt')

%% Make a bar graph of the growth rates for each one

for i = 1:nsubpops
galone(i) = CLLdata(i).g;
gwithk(i) = CLLdata(i).gandK(1);
Ki(i) =CLLdata(i).gandK(2);
end
xvals = 1:1:nsubpops;
figure;
for i = 1:nsubpops
bar(xvals(i),galone(i), 'facecolor', CLLdata(i).color)
hold on
end
set(gca,'XTickLabel',sampsnames)
set(gca,'XTick',8)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('sample')
ylabel('growth rate with set K')
title('Growth rate comparisons with set K')

figure;
for i = 1:nsubpops
bar(xvals(i),gwithk(i), 'facecolor', CLLdata(i).color)
hold on
end
set(gca,'XTick',1:7)
set(gca,'XTickLabel',sampsnames)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('sample')
ylabel('growth rate when fitting K')
title('Growth rate comparisons with individually fit K')
%%
allnames = horzcat(sampsnames, 'global K')

figure;
for i = 1:nsubpops
bar(xvals(i),Ki(i), 'facecolor', CLLdata(i).color)
hold on
bar(8, Kglob)
end
set(gca,'XTick',1:8)
set(gca,'XTickLabel',allnames)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('sample')
ylabel('fit for K')
title('Individually fit K')
%% Consider cutting off the data and setting K and fitting

for i = 1:nsubpops
    
    sigmavec = CLLdata(i).Nstd(1:7);
    sigmavec(1) = 1e4;
    sigma = 1e4; % add a little uncertainty to initial condition
    [p, modelfitg] = fit_logistic(CLLdata(i).Nmean(1:7), tdata(1:7), sigma, Kglob);
    CLLdata(i).gshort = p(1);
    CLLdata(i).modelfitshort = modelfitg;
end
figure;
for i = 1:nsubpops
    subplot(1, nsubpops, i)
      errorbar(CLLdata(i).time(1:7), CLLdata(i).Nmean(1:7),1.96*CLLdata(i).Nstd(1:7), '*', 'color', CLLdata(i).color, 'LineWidth', 2)
         hold on
         plot(CLLdata(i).time(1:7), CLLdata(i).modelfitshort, '-', 'color', CLLdata(i).color, 'LineWidth', 2) 
         %text(CLLdata(i).time(end-2), CLLdata(i).modelfit(end-2), ['g=', num2str(CLLdata(i).gfit)])
         %plot(CLLdata(i).time, CLLdata(i).Nmean + 1.96*CLLdata(i).Nstd, 'color', CLLdata(i).color)
         %plot(CLLdata(i).time, CLLdata(i).Nmean - 1.96*CLLdata(i).Nstd, 'color', CLLdata(i).color)
        xlabel('time (hours)')
        ylabel('N(t)')
        title(['g=', num2str(CLLdata(i).gshort)])
        set(gca,'FontSize',20,'LineWidth',1.5)
        xlim([tdata(1) tdata(7)])
        legend(sampsnames(i), 'Location', 'NorthWest')
        legend boxoff
        ylim([0 6e7])
        
end
%%  Plot fit to first few time points
for i = 1:nsubpops
gshorti(i) = CLLdata(i).gshort;
end
xvals = 1:1:nsubpops;
figure;
for i = 1:nsubpops
bar(xvals(i),gshorti(i), 'facecolor', CLLdata(i).color)
hold on
end
set(gca,'XTick',1:8)
set(gca,'XTickLabel',sampsnames)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('sample')
ylabel('growth rate from shortened data set')
title('Growth rate with shortened data')
%% Fit each inidividual trajectory from shortened time points
for i = 1:nsubpops
    sigmavec(1) = 1e4;
    sigma = 1e4; % add a little uncertainty to initial condition
    for j = 1:3
    [gi, modelfiti] = fit_logistic(CLLdata(i).rawN(1:7, j), tdata(1:7), sigma, Kglob);
    gall(j)= gi;
    modelfits(:,j)= modelfiti;
    end
    CLLdata(i).gall = gall;
    CLLdata(i).modelfitall = modelfits;
end
figure;
for i = 1:nsubpops
    subplot(1, nsubpops, i)
    for j = 1:3
      plot(CLLdata(i).time(1:7), CLLdata(i).rawN(1:7,j), '*', 'color', CLLdata(i).color, 'LineWidth', 2)
         hold on
        plot(CLLdata(i).time(1:7), CLLdata(i).modelfitall(:,j), '-', 'color', CLLdata(i).color, 'LineWidth', 2) 
    end 
        %text(CLLdata(i).time(end-2), CLLdata(i).modelfit(end-2), ['g=', num2str(CLLdata(i).gfit)])
         %plot(CLLdata(i).time, CLLdata(i).Nmean + 1.96*CLLdata(i).Nstd, 'color', CLLdata(i).color)
         %plot(CLLdata(i).time, CLLdata(i).Nmean - 1.96*CLLdata(i).Nstd, 'color', CLLdata(i).color)
        xlabel('time (hours)')
        ylabel('N(t)')
        %title(['g=', num2str(CLLdata(i).gshort)])
        set(gca,'FontSize',20,'LineWidth',1.5)
        xlim([tdata(1) tdata(7)])
        legend(sampsnames(i), 'Location', 'NorthWest')
        legend boxoff
        ylim([0 6e7])
        
end
%%
figure;
for i = 1:nsubpops
    noise =0.1*(0.5-rand(3,1));
    plot((noise+i), CLLdata(i).gall, '*', 'color', CLLdata(i).color, 'LineWidth', 4)
    hold on
end
set(gca,'XTick',1:7)
set(gca,'XTickLabel',sampsnames)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('sample')
ylabel('growth rate')
title('Individual growth rate estimates ')
%% Save the data structure
% this saves the fitted data structure, obtained from the raw data
% structure (run load_raw_data.m)
save('../out/CLLdatagrowth.mat', 'CLLdata')