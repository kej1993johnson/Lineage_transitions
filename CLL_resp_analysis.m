% Response + Relapse of isolated subpopulations of CLL cells

% This script loads in the CLL response and relapse data of the following 
% isolated populations:

% TPO (8 days post thaw, t=0 for scRNAseq)
% TP0 CD18+ (cells from this TP0 population sorted into by high CD18 
% expression levels)
% TP0 CXCR4 (cells from TP0 sorted by high CXCR4 levels (thought to be high
% tolerance)
% FM1 (cells post fluderamine treatment from sample 1)
% FM7 (cells post fluderamine treatment frm sample 7)
% FM1 CD18 (cells post fluderamine treatment sorted for high CD18)
% FM1 CXCR4 (cells from fluderamine treatment sorted for high CXCR4)

close all; clear all; clc;
[resp_data, T]= xlsread('../data/CLL_resp_data.xlsx');
%% Adjust the data
tdays = resp_data(:,1); % this is in days
tdata = tdays.*24; % make it in hours
% cell numbers should be in millions
N_cells = resp_data(:,2:end).*1e6;
nsubpops = 7;
%% Make a data structure that holds each of samples
ntps = length(tdata);
CLLdata = struct('time', zeros(ntps,1), 'rawN', zeros(ntps,3), 'sample', strings(1,1));
    
% Go through and add each sample
sampsnames = {'TP0', 'TP0 CD18+', 'TP0 CXCR4+', 'FM1', 'FM1 CD18+', 'FM1 CXR4+', 'FM7'};
colorsets = varycolor(nsubpops); % give each sample a unique color (for plotting)
for i = 1:nsubpops
    CLLdata(i).time = tdata;
    CLLdata(i).sample = sampsnames(i);
    CLLdata(i).rawN = N_cells(:, i);
    CLLdata(i).color =colorsets(i,:);
end

%% Make a plot of the mean and stdev of data
% Plot the mean of the data
figure;
for i = 1:nsubpops
    plot(CLLdata(i).time, CLLdata(i).rawN,'*-', 'color', CLLdata(i).color, 'LineWidth', 2)
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

sigma = 1e4; % assume our average uncertainty is about 1000 cells?
for j = 1:nsubpops
    ydata = CLLdata(j).rawN;
    N0 = ydata(1);
%     ydata = ydata(2:end);
    ytime = CLLdata(j).time;
%     ytime = ytime(2:end);
    
    % Set up forward models, fit all three nested versions of model
    modelfungood2 = @(p)simmodel2(p,ytime, N0); % double exponential model function with ytime
     
    % INITIAL GUESSES BASED ON DATA
        gguess = log(ydata(end)./ydata(end-1))./(ytime(end)-ytime(end-1))
        phiguess = (ydata(end)/ydata(1)).*exp(-gguess*ytime(end))
        
        if phiguess>0.2
            phiguess = 0.01;
        end
        
        kguess = 0.02;
   
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
    CLLdata(j).params = pbxform2(phatbest2);
    biexpmod = simmodel2(pbxform2(phatbest2), ytime, N0); % mode'
    CLLdata(j).biexpmodel = biexpmod; % add N0
    residuals2 = biexpmod-ydata;
    ybar = mean(ydata);
    Rsq2 = 1-(sum((residuals2).^2)./(sum((ybar-ydata).^2)));
    CLLdata(j).Rsq2 = Rsq2;
    num_params2 = 3;
    n = length(ytime);
    
        ext = 240;
        text = 0:1:ytime(end)+ext;
        CLLdata(j).biexpextend = simmodel2(pbxform2(phatbest2), text, N0); % model
        CLLdata(j).tmod = text;
        icrit = find(CLLdata(j).biexpextend>1.2*N0, 1, 'first');
        CLLdata(j).TTB120= text(icrit);
 end
 %% Flip through model fits
 
figure;

for j = 1:nsubpops-1
    subplot(2,3,j)
    plot(CLLdata(j).time, CLLdata(j).rawN, '*-', 'LineWidth', 2, 'color', num2str(CLLdata(j).color))
    hold on
    plot(CLLdata(j).tmod, CLLdata(j).biexpextend, 'LineWidth', 2, 'color',num2str(CLLdata(j).color))   
    plot(CLLdata(j).TTB120, 1.2*N0, 'k*', 'LineWidth', 5)
    plot([CLLdata(j).TTB120, CLLdata(j).TTB120], [0, 1.5*N0], 'k--', 'LineWidth', 4)
    %text(CLLdata(j).TTB120, 1.2*N0, 'TTB120')
    title([CLLdata(j).sample])
    set(gca,'FontSize',20,'LineWidth',1.5)
    legend('data', 'model', ['TTB120= ', num2str(CLLdata(j).TTB120),' days'], 'FontSize', 14, 'Location', 'NorthWest')
    legend boxoff
    xlabel ('time (hours)')
    ylabel('N(t)')
    xlim([0 CLLdata(j).TTB120+24])
    ylim([ 0 1.5*N0])
end

%% Save the loaded data 
save('../out/CLLdataresp.mat', 'CLLdata')
%% Make bar graph of TTB120 
TTB120s = [];
phis = [];
ks = [];
gs = [];
for i = 1:nsubpops
TTB120s=vertcat(TTB120s, CLLdata(i).TTB120);
phis = vertcat(phis, CLLdata(i).params(1));
gs = vertcat(gs, CLLdata(i).params(2));
ks = vertcat(ks, CLLdata(i).params(3));
end
xvals = 1:1:nsubpops;
figure;
for i = 1:nsubpops
bar(xvals(i),TTB120s(i), 'facecolor', CLLdata(i).color)
hold on
end
set(gca,'XTickLabel',sampsnames)
set(gca,'XTick',1:7)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('sample')
ylabel('TTB120')
title('Time to 1.2*N_{0}')
%% Plot things
% Plot resistant fraction
figure;
for i = 1:nsubpops
bar(xvals(i),phis(i), 'facecolor', CLLdata(i).color)
hold on
end
set(gca,'XTickLabel',sampsnames)
set(gca,'XTick',1:7)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('sample')
ylabel('Fraction resistant')
title('Fraction of cells resistant to drug at onset')
% Plot deathr ate due to drug (k)
figure;
for i = 1:nsubpops
bar(xvals(i),ks(i), 'facecolor', CLLdata(i).color)
hold on
end
set(gca,'XTickLabel',sampsnames)
set(gca,'XTick',1:7)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('sample')
ylabel('death rate due to drug (k)')
title('Death rate due to drug (k)')
% Plot resistant regrowth rate
figure;
for i = 1:nsubpops
bar(xvals(i),gs(i), 'facecolor', CLLdata(i).color)
hold on
end
set(gca,'XTickLabel',sampsnames)
set(gca,'XTick',1:7)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('sample')
ylabel('resistant cell growth rate')
title('Resistant cell regrowth rate')
