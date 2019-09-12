% Barcode sampling 
% This script loads in the barcode sampling data which contains raw numbers
% of barcodes found in samples at three different time points after initial
% barcoding:
%OGTP1 = 48h, 3% population which was 19.8x10^6 cells
%OGTP2 = 96h, 3% population which was 75.9x10^6 cells 
%OGTP3 = 144h, 3% population was at 180.8x10^6 cells


close all; clear all; clc
%% Load in OG and TPO barcode frequency distribution
[N, T]= xlsread('../data/Barcode_sampling_pre_treat.xlsx');
% Want to make a structure that tracks the number of barcodes, the barcode
% the rank based on the first time point, second time point, and last time
% point and the scaled population number

% Want to set all NaNs to 0 in the N matrix
idnotobs = isnan(N);
N(isnan(N))=0;

%% Make a structure that holds the barcode information
% Initialize the structure
% For each barcode have a vector of times, a vector of number of cells at
% that time whos barcode was captured, and then the number of cells that
% corresponds to if we extrapolate to the whole expanding population, and
% then the proportion of cells in the population at that time
nCD18 = 431200; 
nCXCR4 = 71867;
numsamps = [0.594e6, 2.277e6, 5.42e6, nCD18, nCXCR4]; % number of cells loaded in each sample submitted for barcoding
numtots = [ 19.8e6, 75.9e6, 180.8e6, 10*nCD18, 10*nCXCR4]; % total number of cells
nreads = sum(N); % sum each column
pct_samps = numsamps(1:3)./numtots(1:3); % these should all be around 3 but just to check
% To estimate the number of cells in each barcode in the sample and the
% population you need
bcdabund = N./nreads; % proportion of barcodes at each time point
Nsamps = bcdabund.*numsamps; % number of cells with that barcode in the sample
Ntot = bcdabund.*numtots; % number of cells with that barcode at each time point
%%
bcd = struct( 'time', zeros(1,1),'barcode',strings(1), 'nsamp',...
    zeros(1,1),'ntot', zeros(1,1),'abund', zeros(1,1), 'rank1',...
    zeros(1,1), 'rank2', zeros(1,1), 'rank3', zeros(1,1));
% Make the rankt1, rankt2, and rankt3 list of indices
% Use these to keep track of how the barcode moves in ranking
[B1, rankt1]=sort(N(:,1), 'descend');
[B2, rankt2] = sort(N(:,2), 'descend');
[B3, rankt3] = sort(N(:,3), 'descend');
%% Loop through and add each value to the barcode structure
num_bcds = length(N);
for i = 1:num_bcds
    bcd(i).time = [48, 96, 144];
    bcd(i).barcode = char(T(i+1,1));
    bcd(i).nsamp = Nsamps(i,:);
    bcd(i).ntot = Ntot(i,:);
    bcd(i).abund = bcdabund(i,:);
    bcd(i).rank1 = rankt1(i);
    bcd(i).rank2 =rankt2(i);
    bcd(i).rank3 = rankt3(i);
    % flag the barcode if it is not fully observed in all three TPs
    if any(idnotobs(i,:))
    bcd(i).allobs = 0;
    else
        bcd(i).allobs = 1;
    end
end
%% Look at the barcode distribution histogram
figure;
subplot(3,1,1)
bar(1:1:num_bcds, (bcdabund(:,1)))
xlabel('barcode')
ylabel('abundance')
title('t=0h barcode abundance distribution')
ylim([0 0.1])
%xlim([0, 1000])
set(gca,'FontSize',20,'LineWidth',1.5, 'YScale', 'log')

subplot(3,1,2)
bar(1:1:num_bcds, (bcdabund(:,2)))
xlabel('barcode')
ylabel('abundance')
title('t=48h ordered by t=0h')
ylim([0 0.1])
%xlim([0, 1000])
set(gca,'FontSize',20,'LineWidth',1.5, 'YScale', 'log')

subplot(3,1,3)
bar(1:1:num_bcds, (bcdabund(:,3)))
xlabel('barcode')
ylabel('abundance')
title('t=96h ordered by t=0h')
ylim([0 0.1])
%xlim([0, 1000])
set(gca,'FontSize',20,'LineWidth',1.5, 'YScale', 'log')

% Make a heatmap
% make this into a 3 column heat map with barcode on the y-axis and time
% point on the x axis, colored by bcdabundance.
xlabs = {'t=0h', 't=48h', 't=96h'};
figure;
colormap('hot')
imagesc(log10(bcdabund(:, 1:3)))
h = colorbar;
h.Label.String = 'Log10 barcode abundance'
ylabel(h, 'Log10 barcode abundance')
set(gca,'FontSize',16,'LineWidth',1.5)
ylabel('barcode')
title('Heatmap of barcode abundance')
set(gca,'XTick',[1, 2 ,3], 'XTickLabel', xlabs);% something like this
%%
figure;
subplot(3,1,1)
bar(1:1:num_bcds, bcdabund(:,1))
xlabel('barcode')
ylabel('abundance')
title('t=48h barcode abundance distribution')
set(gca,'FontSize',20,'LineWidth',1.5, 'YScale', 'log')
%xlim([0, 1000])
ylim([0 0.097])

subplot(3,1,2)
bar(1:1:num_bcds, bcdabund(rankt2,2))
xlabel('barcode')
ylabel('abundance')
title('t=96h descending order')
set(gca,'FontSize',20,'LineWidth',1.5, 'YScale', 'log')
%xlim([0, 1000])
ylim([0 0.097])

subplot(3,1,3)
bar(1:1:num_bcds, bcdabund(rankt3,3))
xlabel('barcode')
ylabel('abundance')
title('t=144h descending order')
set(gca,'FontSize',20,'LineWidth',1.5,'YScale', 'log')
%xlim([0, 1000])
ylim([0 0.097])
%% Plot as separate figures
figure;
bar(1:1:num_bcds, bcdabund(:,1))
xlabel('barcode')
ylabel('abundance')
title('t=48h')
set(gca,'FontSize',20,'LineWidth',1.5,'YScale', 'log')
%xlim([0, 500])
ylim([0 0.097])

figure;
bar(1:1:num_bcds, bcdabund(rankt2,2))
xlabel('barcode')
ylabel('abundance')
title('t=96h ')
set(gca,'FontSize',20,'LineWidth',1.5, 'YScale', 'log')
%xlim([0, 500])
ylim([0 0.097])

figure;
bar(1:1:num_bcds, bcdabund(rankt3,3))
xlabel('barcode')
ylabel('abundance')
title('t=144h')
set(gca,'FontSize',20,'LineWidth',1.5, 'YScale', 'log')
%xlim([0, 500])
ylim([0 0.097])

%% Plot the indiviudal growth trajectories
figure;
for i = 1:1000%num_bcds
    tvec =  bcd(i).time;
    nvec =  bcd(i).ntot;
    %tvec =  bcd(i).time;
    %nvec =  bcd(i).ntot;
    plot(tvec, log(nvec/nvec(1)), '-', 'LineWidth', 1)
    set(gca,'FontSize',20,'LineWidth',1.5)
    hold on
    xlabel('time (hours)')
    ylabel('log(N(t)/N_{0})')
    title('Normalized growth rate of each lineage')
    xlim([tvec(1) tvec(end)])
end

figure;
for i = 1:1000%num_bcds
    tvec =  bcd(i).time;
    nvec =  bcd(i).ntot;
    %tvec =  bcd(i).time;
    %nvec =  bcd(i).ntot;
    plot(tvec, nvec, '-', 'LineWidth', 1)
    set(gca,'FontSize',20,'LineWidth',1.5, 'Yscale', 'log')
    hold on
    xlabel('time (hours)')
    ylabel('N(t)')
    title('Raw lineage growth trajectories in log scale')
    xlim([tvec(1) tvec(end)])
end

figure;
for i = 1:1000%num_bcds
    tvec =  bcd(i).time;
    nvec =  bcd(i).ntot;
    %tvec =  bcd(i).time;
    %nvec =  bcd(i).ntot;
    plot(tvec, nvec, '-', 'LineWidth', 1)
    set(gca,'FontSize',20,'LineWidth',1.5)
    hold on
    xlabel('time (hours)')
    ylabel('N(t)')
    title('Raw lineage growth trajectories')
    xlim([tvec(1) tvec(end)])
end
%% Fit the trajectories from the 3 observed time points for each lineage



for i = 1:num_bcds
  if bcd(i).allobs ==1
      tvec = bcd(i).time - bcd(i).time(1); % shift time points for fitting set t=48h to t=0 for now
      nvec = bcd(i).ntot(2:end);
      num_params = 1;
      % Fit on all but t=48h, which we will call t=0
      N0 = bcd(i).ntot(1); % set N0 as the first observed tp at t=48h
      nbar = mean(nvec); % average value of volumes over all measured times
      
      % First fit each curve to single exponential over entire time course
      LB = -Inf ;  % Lower Bounds
      UB = Inf; % Upper Bounds
      params0 = log(nvec(2)/nvec(1))./(tvec(2)-tvec(1)); % best guess (slope of the last two tps)
      options = optimset('TolFun',1e-12,'Tolx',1e-12,'MaxIter',1000,'Display','off','FinDiffRelStep',1e-3);
      [g, resnorm, reslsq]= lsqnonlin(@fitsingleexp, params0, LB, UB, options, nvec, N0, tvec);
      bcd(i).Rsq = 1- (sum((reslsq.^2))./(sum((nbar-nvec).^2)));
      bcd(i).Nmodel = singleexpmodel(g,N0, tvec); % model fit for plotting (include t=0)
      % Only keep g if it's a decent fit
      if bcd(i).Rsq >0.4
          bcd(i).g = g;
          bcd(i).badfit = 0;
      else % keep track of bad fits
          bcd(i).badfit = 1;
          bcd(i).g = g;
      end
  else % if not all three time points are observed, don't try to fit
      bcd(i).Rsq = 0;
      bcd(i).Nmodel = zeros(length(bcd(i).time), 1);
      bcd(i).g = 0;
      bcd(i).badfit = 1;
  end
end
%% Plot some example results of fit versus data
figure;
for i = 1:10
    plot(bcd(i).time, bcd(i).ntot, '*')
    hold on
    plot(bcd(i).time, bcd(i).Nmodel, '-', 'LineWidth', 2)
    set(gca,'FontSize',20,'LineWidth',1.5)
    xlabel('time')
    ylabel('N(t)')
    legend('data','fit', 'Location', 'Northwest')
    legend boxoff
    title(['g=', num2str(round(bcd(i).g,4))])
    pause
end
%% Plot growth rate distribution
glist = [];
for i = 1:num_bcds
    if bcd(i).allobs ==1 && bcd(i).badfit ==0
    glist = vertcat(glist, bcd(i).g);
    end
end
[gdistrib, rankg] = sort(glist, 'descend');
gmean = mean(glist)
gmed = median(glist)
gstdev = std(glist)
%% Plot growth rate distribution
figure;
bar(1:1:length(glist), glist)
xlabel('barcode')
ylabel('growth rate')
title('growth rate distribution')
set(gca,'FontSize',20,'LineWidth',1.5)
%xlim([0, 500])
%ylim([0 0.097])

figure;
bar(1:1:length(glist), gdistrib)
xlabel('barcode')
ylabel('growth rate')
title('growth rate distribution')
set(gca,'FontSize',20,'LineWidth',1.5)
%xlim([0, 500])
%ylim([0 0.097])

figure
histogram(glist)
xlabel('growth rate')
ylabel('number of lineages')
title('growth rate distribution')
set(gca,'FontSize',20,'LineWidth',1.5)