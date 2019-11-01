% Response + Relapse of isolated subpopulations of CLL cells

% This script loads in the CLL response and relapse data of the following 
% isolated populations in triplicate

% TPO (8 days post thaw, t=0 for scRNAseq)
% TP0 CD18+ (cells from this TP0 population sorted into by high CD18 
% expression levels)
% TP0 CXCR4 (cells from TP0 sorted by high CXCR4 levels (thought to be high
% tolerance)

% Only the TP0 CD18 psoitive are monitored for the full 26 days


close all; clear all; clc;
[resp_data, T]= xlsread('../data/CLL_resp_relapse_rep2.xlsx');
%% Adjust the data
tdays = resp_data(:,1); % this is in days
tdata = tdays.*24; % make it in hours
% cell numbers should be in millions
N_cells = resp_data(:,2:end).*1e6;
nsubpops = 3;
%% Make a data structure that holds each of samples
ntps = length(tdata);
CLLdata = struct('time', zeros(ntps,1), 'rawN', zeros(ntps,3), 'sample', strings(1,1));
    
% Go through and add each sample
sampsnames = {'TP0', 'TP0 CD18+', 'TP0 CXCR4+'};
colorsets = varycolor(nsubpops); % give each sample a unique color (for plotting)

% load in the TP0 data sets
for i = 1:3
    CLLdata(i).time = 24*tdata(1:end-4,1);
    CLLdata(i).sample = sampsnames(1);
    CLLdata(i).rawN = N_cells(1:end-4, i);
    CLLdata(i).color =colorsets(1,:);
end

% load in the TP0CD18 data set
for i = 1:3
    CLLdata(i+3).time = 24*tdata(1:end-3);
    CLLdata(i+3).sample = sampsnames(2);
    CLLdata(i+3).rawN = N_cells(1:end-3, i+3);
    CLLdata(i+3).color =colorsets(2,:);
end

% load in the TP0 CXCR4 data set
for i = 1:3
    CLLdata(i+6).time = 24*tdata(1:end-5);
    CLLdata(i+6).sample = sampsnames(3);
    CLLdata(i+6).rawN = N_cells(1:end-5, i+6);
    CLLdata(i+6).color =colorsets(3,:);
end
%% Make a plot of the raw data
% Plot the mean of the data
figure;
for i = 1:length(CLLdata)
    plot(CLLdata(i).time, CLLdata(i).rawN,'*-', 'color', CLLdata(i).color, 'LineWidth', 2)
         hold on
         %text(CLLdata(i).time(end-2), CLLdata(i).Nmean(end-2), [(CLLdata(i).sample)])
         %plot(CLLdata(i).time, CLLdata(i).Nmean + 1.96*CLLdata(i).Nstd, 'color', CLLdata(i).color)
         %plot(CLLdata(i).time, CLLdata(i).Nmean - 1.96*CLLdata(i).Nstd, 'color', CLLdata(i).color)
        xlabel('time (hours)')
        ylabel('N(t)')
        title('N(t) for each sample')
        set(gca,'FontSize',20,'LineWidth',1.5)
        %xlim([tdata(1) tdata(end)])
       
end
%% Make CLL sum structure

for i = 1:nsubpops
    CLLsum(i).time = [];
    CLLsum(i). Nmat = [];
    CLLsum(i).Nmean = [];
    CLLsum(i).Nstd = [];
    CLLsum(i).color = colorsets(i,:);
    CLLsum(i).sample = sampsnames(i);
end
% load in TP0
CLLsum(1).time = CLLdata(1).time;
for i = 1:3
CLLsum(1).Nmat = horzcat(CLLsum(1).Nmat, CLLdata(i).rawN);
end
CLLsum(1).Nmean = mean(CLLsum(1).Nmat,2);
CLLsum(1).Nstd = std(CLLsum(1).Nmat,0, 2);

% load in TP0 CD18
CLLsum(2).time = CLLdata(4).time;
for i = 4:6
CLLsum(2).Nmat = horzcat(CLLsum(2).Nmat, CLLdata(i).rawN);
end
CLLsum(2).Nmean = mean(CLLsum(2).Nmat,2);
CLLsum(2).Nstd = std(CLLsum(2).Nmat,0, 2);

% load in TP0 CXCR4
CLLsum(3).time = CLLdata(7).time;
for i = 7:9
CLLsum(3).Nmat = horzcat(CLLsum(3).Nmat, CLLdata(i).rawN);
end
CLLsum(3).Nmean = mean(CLLsum(3).Nmat,2);
CLLsum(3).Nstd = std(CLLsum(3).Nmat,0, 2);
% Plot the mean of the data
figure;
for i = 1:length(CLLsum)
        errorbar(CLLsum(i).time, CLLsum(i).Nmean,0.5.*1.96.*CLLsum(i).Nstd, '*-', 'color', CLLsum(i).color, 'LineWidth', 2)
         hold on
         %text(CLLdata(i).time(end-2), CLLdata(i).Nmean(end-2), [(CLLdata(i).sample)])
         %plot(CLLdata(i).time, CLLdata(i).Nmean + 1.96*CLLdata(i).Nstd, 'color', CLLdata(i).color)
         %plot(CLLdata(i).time, CLLdata(i).Nmean - 1.96*CLLdata(i).Nstd, 'color', CLLdata(i).color)
         legend ('TP0', 'TP0 CD18', 'TP0 CXCR4', 'Location', 'NorthWest')
         legend boxoff
        xlabel('time (hours)')
        ylabel('N(t)')
        title('N(t) for each sample')
        set(gca,'FontSize',20,'LineWidth',1.5)
        %xlim([tdata(1) tdata(end)])
       
end
%% Fit each trajectory to the bi-exponential model of response and relapse
% The model looks like this
%
% $$ N(t) = N_0 [ \phi e^{gt} + (1-\phi) e^{-kt} ] $$
% 
% where:
%
% * $N_0$ is the initial cell number
% * $\phi$ is the initial resistant fraction
% * g>0 is the resistant growth rate
% * k>0 is the kill rate (actually -net "growth rate" on treatment) on sensitive cells
%

% Define transforms 
yfxform = @(y)log(y); % 'forward' transform for data and model output
ybxform = @(yhat)exp(yhat); % 'inverse' transform for data and model output
% double exponential
pfxform2 = @(pval)[0 1 1].*log(pval)+[1 0 0].*log(pval./(1-pval)); %'forward' parameter transform into Reals
pbxform2 = @(phat)[0 1 1].*exp(phat)+[1 0 0].*(exp(phat)./(1+exp(phat)));  %'backward' parameter transform into model space

 % assume our average uncertainty is about 1000 cells?
for j = 1:nsubpops
    ydata = CLLsum(j).Nmean;
    sigma = CLLsum(j).Nstd;
    N0 = ydata(1);
%     ydata = ydata(2:end);
    ytime = CLLsum(j).time;
%     ytime = ytime(2:end);
    
    % Set up forward models, fit all three nested versions of model
    modelfungood2 = @(p)simmodel2(p,ytime, N0); % double exponential model function with ytime
     
    % INITIAL GUESSES BASED ON DATA
        gguess = log(ydata(end)./ydata(end-2))./(ytime(end)-ytime(end-2))
        phiguess = (ydata(end)/ydata(1)).*exp(-gguess*ytime(end))
        
        if phiguess>0.2
            phiguess = 0.01;
        end
        
        kguess = log(ydata(1)./ydata(2))./ytime(2)-ytime(1);
   
    theta2 = [phiguess, gguess, kguess];
    
    % Write log likelihood function based on assumption of normally
    % distributed sampling error
    
    % Goal: maximize the probability of the data given the model. Normpdf
    % will output a probability of the data (x- 1st argument), given the
    % mean(expectation, here the model), and the variance, at each time point. 
    % take the log and minimize the NLL to maximize likeihood of data given
    % the model
  
    % double exponential
    loglikelihood2 = @(phat)sum(log(normpdf(yfxform(ydata),yfxform(modelfungood2(pbxform2(phat))), sigma)));
  
    % Write objective functions for each model
    objfun2 = @(phat)-loglikelihood2(phat);
    phatbest2 = fminsearch(objfun2, pfxform2(theta2));
    
    % Save some stuff
    % save biexp stuff (most relevant so will put first)
    CLLsum(j).params = pbxform2(phatbest2);
    biexpmod = simmodel2(pbxform2(phatbest2), ytime, N0); % mode'
    CLLsum(j).biexpmodel = biexpmod; % add N0
    residuals2 = biexpmod-ydata;
    ybar = mean(ydata);
    Rsq2 = 1-(sum((residuals2).^2)./(sum((ybar-ydata).^2)));
    CLLsum(j).Rsq2 = Rsq2;
    num_params2 = 3;
    n = length(ytime);
    
        ext = 3000;
        text = 0:1:ytime(end)+ext;
        CLLsum(j).biexpextend = simmodel2(pbxform2(phatbest2), text, N0); % model
        CLLsum(j).tmod = text;
        icrit = find(CLLsum(j).biexpextend>1.2*N0, 1, 'first');
        CLLsum(j).TTB120= text(icrit);
 end
 %% Flip through model fits
 
figure;

for j = 1:nsubpops
    subplot(1,3,j)
    errorbar(CLLsum(j).time, CLLsum(j).Nmean, 0.5*1.96*CLLsum(j).Nstd, '*-', 'LineWidth', 2, 'color', num2str(CLLsum(j).color))
    hold on
    plot(CLLsum(j).tmod, CLLsum(j).biexpextend, 'LineWidth', 2, 'color',num2str(CLLsum(j).color))   
    %plot([0;text(end)], [1.2*N0; 1.2*N0], 'r--', 'LineWidth', 2)
     plot(CLLsum(j).TTB120, 1.2*N0, 'r*', 'LineWidth', 5)
    %plot([CLLdata(j).TTB120, CLLdata(j).TTB120], [0, 1.5*N0], 'k--', 'LineWidth', 4)
    %text(CLLdata(j).TTB120, 1.2*N0, 'TTB120')
    title([CLLsum(j).sample])
    set(gca,'FontSize',20,'LineWidth',1.5)
    legend('data', 'model', ['TTB120= ', num2str(CLLsum(j).TTB120),' hours'], 'FontSize', 14, 'Location', 'NorthWest')
    legend boxoff
    xlabel ('time (hours)')
    ylabel('N(t)')
    %xlim([0 CLLdata(j).TTB120+24])
    %ylim([ 0 1.5*N0])
end

%% Save the loaded data 
save('../out/CLLdataresp2.mat', 'CLLsum')
%% Make bar graph of TTB120 
TTB120s = [];
phis = [];
ks = [];
gs = [];
for i = 1:nsubpops
TTB120s=vertcat(TTB120s, CLLsum(i).TTB120);
phis = vertcat(phis, CLLsum(i).params(1));
gs = vertcat(gs, CLLsum(i).params(2));
ks = vertcat(ks, CLLsum(i).params(3));
end
xvals = 1:1:nsubpops;
figure;
for i = 1:nsubpops
bar(xvals(i),TTB120s(i), 'facecolor', CLLsum(i).color)
hold on
end
set(gca,'XTickLabel',sampsnames)
set(gca,'XTick',1:3)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('sample')
ylabel('TTB120')
title('Time to 1.2*N_{0}')
%% Plot things
% Plot resistant fraction
figure;
for i = 1:nsubpops
bar(xvals(i),phis(i), 'facecolor', CLLsum(i).color)
hold on
end
set(gca,'XTickLabel',sampsnames)
set(gca,'XTick',1:3)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('sample')
ylabel('Fraction resistant')
title('Fraction of cells resistant to drug at onset')
% Plot deathr ate due to drug (k)
figure;
for i = 1:nsubpops
bar(xvals(i),ks(i), 'facecolor', CLLsum(i).color)
hold on
end
set(gca,'XTickLabel',sampsnames)
set(gca,'XTick',1:3)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('sample')
ylabel('death rate due to drug (k)')
title('Death rate due to drug (k)')
% Plot resistant regrowth rate
figure;
for i = 1:nsubpops
bar(xvals(i),gs(i), 'facecolor', CLLsum(i).color)
hold on
end
set(gca,'XTickLabel',sampsnames)
set(gca,'XTick',1:3)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('sample')
ylabel('resistant cell growth rate')
title('Resistant cell regrowth rate')
