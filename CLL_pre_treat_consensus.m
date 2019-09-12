% Barcode sampling 
% This script loads in the barcode sampling data which contains raw numbers
% of barcodes found in samples at three different time points after initial
% barcoding:
%OGTP1 = 48h, 3% population which was 19.8x10^6 cells
%OGTP2 = 96h, 3% population which was 75.9x10^6 cells 
%OGTP3 = 144h, 3% population was at 180.8x10^6 cells

% For the CD18 and CXCR4 populations the number of cells sampled is:
% CD18 cells = 431,200 % sensitive
% CXCR4 cells = 71,867 % resistant

% This is the same script as the CLL_pre_treat_barcode_sampling but now we
% generate the distributions and we remove any barcodes who don't have
% three observations, making a list
close all; clear all; clc
%% Load in OG and TPO barcode frequency distribution
[Np, T]= xlsread('../data/Barcode_sampling_pre_treat.xlsx');
[Npre, Tpre] = xlsread('../data/barcode_sampling_t0hr.xlsx');
%% Find number of reads in each column before removing data 
nreadsp = nansum(Np); % sum each column
nreads0 = nansum(Npre);
nreads = horzcat(nreads0,nreadsp);

Npre(isnan(Npre))=0;
%% Combine the t=0 hr data set with the t=48hr data set
Tpost = cell2table(T(2:end,1));
Tpost.t48 = Np(:,1);
Tpost.t96 = Np(:,2);
Tpost.t144 = Np(:,3);
Tpost.CD18pos = Np(:,4);
Tpost.CXCR4post = Np(:,5);
Tzero = cell2table(Tpre(2:end,1));
Tzero.t0 = Npre;
%% Make your consensus and all barcode tables
% Make a table that joins based on the Tzero barcodes
[Tcons, ia, ib] = outerjoin(Tzero, Tpost, 'Type', 'left', 'MergeKeys', true); 
% This should join by the key variable in common, which is Var1.

% Make a table that adds all of the barcodes (unique from both Tzero and
% Tpost
Tall = outerjoin(Tzero, Tpost, 'MergeKeys', true); 
% This should join those barcodes that are the same and then combine ones that are different
%% Find the indexes from the table that order by t=0hr
N0 = Tcons.t0;
%N0(isnan(N0))=0;
[check,iord]=sort(N0, 'descend');
Tsort = sortrows(Tcons,{'t0'}, 'descend');

% Set NaNs to 0 in the first time point
Tall.t0(isnan(Tall.t0))=0;
% Once NaNs =0, sort by T0
Tsortall = sortrows(Tall,{'t0'}, 'descend');



%% Reset your Ns
Ntbl = Tsort(:,2:end);
Ncell = table2cell(Ntbl);
N = cell2mat(Ncell); % This one is 13648 unique barcodes

Ntblall = Tsortall(:,2:end);
Ncellall = table2cell(Ntblall);
Nall = cell2mat(Ncellall);
Nall(isnan(Nall))=0; % this one is 17730 unique barcodes


%% Remove barcodes not observed in all 4 time points
% Want to make a structure that tracks the number of barcodes, the barcode
% the rank based on the first time point, second time point, and last time
% point and the scaled population number

% Want to set all NaNs to 0 in the N matrix
idnotobs = isnan(N(:,1:4)); % We want only barcodes observed in all four time points 
N(isnan(N))=0;
Ncons=[];
cons_barcodes = {};
for i = 1:length(N)
    if ~any(idnotobs(i,:))
        Ncons = vertcat(Ncons, N(i,:));
        cons_barcodes = vertcat(cons_barcodes, char(Tsort.Var1(i)));
    end
end


%% Make a structure that holds the consensus barcode information
% Initialize the structure
% For each barcode have a vector of times, a vector of number of cells at
% that time whos barcode was captured, and then the number of cells that
% corresponds to if we extrapolate to the whole expanding population, and
% then the proportion of cells in the population at that time
n0temps = 1e6; % THIS IS WRONG
nCD18 = 431200; 
nCXCR4 = 71867;
numsamps = [n0temps, 0.594e6, 2.277e6, 5.42e6, nCD18, nCXCR4]; % number of cells loaded in each sample submitted for barcoding
numtots = [10e6, 19.8e6, 75.9e6, 180.8e6, 10*nCD18, 10*nCXCR4]; % total number of cells 
pct_samps = numsamps(2:4)./numtots(2:4); % these should all be around 3 but just to check
% To estimate the number of cells in each barcode in the sample and the
% population you need
% consensus 7532 barcodes
bcdabund = Ncons./nreads; % proportion of barcodes at each time point
Nsamps = bcdabund.*numsamps; % number of cells with that barcode in the sample
Ntot = bcdabund.*numtots; % number of cells with that barcode at each time point

% all barcodes 17730 barcodes
bcdabundall = Nall./nreads;
Nsampsall = bcdabundall.*numsamps;
Ntotall = bcdabundall.*numtots;
%% Set up consensus structure
bcd = struct( 'time', zeros(1,1),'barcode',strings(1), 'nsamp',...
    zeros(1,1),'ntot', zeros(1,1),'abund', zeros(1,1),'CD18pos', zeros(1,1), 'CXCR4pos', zeros(1,1),'rank1',...
    zeros(1,1), 'rank2', zeros(1,1), 'rank3', zeros(1,1), 'rank4', zeros(1,1));
% Make the rankt1, rankt2, and rankt3 list of indices
% Use these to keep track of how the barcode moves in ranking
[B1, rankt1]=sort(Ncons(:,1), 'descend');
[B2, rankt2] = sort(Ncons(:,2), 'descend');
[B3, rankt3] = sort(Ncons(:,3), 'descend');
[B4, rankt4] = sort(Ncons(:,4), 'descend');
%% Loop through and add each value to the consensus barcode structure
num_bcds = length(Ncons);
thres_ratio = 1; % abundance to call a lineage CD18 positive
CD18ratio = [];
ct_CXCR4pos = 0;
for i = 1:num_bcds
    bcd(i).time = [0, 48, 96, 144];
    bcd(i).barcode = char(cons_barcodes(i));
    bcd(i).nsamp = Nsamps(i,:);
    bcd(i).ntot = Ntot(i,:);
    bcd(i).abund = bcdabund(i,:);
    bcd(i).rank1 = rankt1(i);
    bcd(i).rank2 =rankt2(i);
    bcd(i).rank3 = rankt3(i);
    bcd(i).rank4 = rankt4(i);
    % flag the barcode if it is not fully observed in all three TPs
    bcd(i).CD18ratio = bcd(i).abund(1,5)/bcd(i).abund(1,6);
    CD18ratio = vertcat(CD18ratio, bcd(i).CD18ratio);
    if bcd(i).CD18ratio >1
        bcd(i).CD18pos = 1;
        bcd(i).CXCR4pos = 0;
    end
    if bcd(i).CD18ratio <1
        bcd(i).CD18pos = 0;
        bcd(i).CXCR4pos = 1;
        ct_CXCR4pos= ct_CXCR4pos +1
    end
    if isnan(bcd(i).CD18ratio)
        bcd(i).CD18pos = 0;
        bcd(i).CXCR4pos = 0;
    end
%     if bcd(i).CD18ratio>thres_ratio
%         bcd(i).CD18pos = 1;
%     else
%         bcd(i).CD18pos = 0;
%     end
%     
%     if bcd(i).abund(1,6)>thres_abund
%         bcd(i).CXCR4pos = 1;
%     else
%         bcd(i).CXCR4pos = 0;
%     end
%     
%     if bcd(i).CD18pos==1 && bcd(i).CXCR4pos ==1
%         if bcd(i).abund(1,5)> bcd(i).abund(1,6) % indicates more of the lineage is in CD18pos
%             bcd(i).CXCR4pos = 0;
%         end
%         if bcd(i).abund(1,5)<bcd(i).abund(1,6) % indicates more of the lineage is in CXCR4pos
%             bcd(i).CD18pos = 0;
%         end
%     end
%     
    
end

pct_CXCR4pos = ct_CXCR4pos/num_bcds
%% Plot the CD18 ratio

CD18ratio(isinf(CD18ratio))=1000;
figure;
histogram(CD18ratio, 1000)
xlabel('CD18/CXCR4 abundance ratio')
ylabel('frequency')
title('CD18/CXCR4 distribution of ratios')
%xlim([0,2])
set(gca,'FontSize',20,'LineWidth',1.5, 'Yscale', 'log')
%% Set up all barcodes structure
bcdall = struct( 'time', zeros(1,1),'barcode',strings(1), 'nsamp',...
    zeros(1,1),'ntot', zeros(1,1),'abund', zeros(1,1),'rank1',...
    zeros(1,1), 'rank2', zeros(1,1), 'rank3', zeros(1,1), 'rank4', zeros(1,1));
% Make the rankt1, rankt2, and rankt3 list of indices
% Use these to keep track of how the barcode moves in ranking
[B1all, rankt1all]=sort(Nall(:,1), 'descend');
[B2all, rankt2all] = sort(Nall(:,2), 'descend');
[B3all, rankt3all] = sort(Nall(:,3), 'descend');
[B4all, rankt4all] = sort(Nall(:,4), 'descend');
%% Loop through and add each value to the ALL barcode structure
num_bcds_all = length(Nall);
for i = 1:num_bcds_all
    bcdall(i).time = [0, 48, 96, 144];
    bcdall(i).barcode = char(Tall.Var1(i));
    bcdall(i).nsamp = Nsampsall(i,:);
    bcdall(i).ntot = Ntotall(i,:);
    bcdall(i).abund = bcdabundall(i,:);
    bcdall(i).rank1 = rankt1all(i);
    bcdall(i).rank2 =rankt2all(i);
    bcdall(i).rank3 = rankt3all(i);
    bcdall(i).rank4 = rankt4all(i);
    % flag the barcode if it is not fully observed in all three TPs
    
end
%% Look at the barcode distribution histogram
figure;
subplot(4,1,1)
bar(1:1:num_bcds, (bcdabund(:,1)))
xlabel('barcode')
ylabel('abundance')
title('t=0h barcode abundance distribution')
ylim([0 0.097])
%xlim([0, 1000])
set(gca,'FontSize',20,'LineWidth',1.5, 'YScale', 'log')

subplot(4,1,2)
bar(1:1:num_bcds, (bcdabund(:,2)))
xlabel('barcode')
ylabel('abundance')
title('t=48h ordered by t=0h')
ylim([0 0.097])
%xlim([0, 1000])
set(gca,'FontSize',20,'LineWidth',1.5, 'YScale', 'log')

subplot(4,1,3)
bar(1:1:num_bcds, (bcdabund(:,3)))
xlabel('barcode')
ylabel('abundance')
title('t=96h ordered by t=0h')
ylim([0 0.097])
%xlim([0, 1000])
set(gca,'FontSize',20,'LineWidth',1.5, 'YScale', 'log')

subplot(4,1,4)
bar(1:1:num_bcds, (bcdabund(:,4)))
xlabel('barcode')
ylabel('abundance')
title('t=144h ordered by t=0h')
ylim([0 0.097])
%xlim([0, 1000])
set(gca,'FontSize',20,'LineWidth',1.5, 'YScale', 'log')
%%
figure;
subplot(4,1,1)
bar(1:1:num_bcds, bcdabund(:,1))
xlabel('barcode')
ylabel('abundance')
title('t=0h barcode abundance distribution')
%set(gca,'FontSize',20,'LineWidth',1.5)
set(gca,'FontSize',20,'LineWidth',1.5, 'YScale', 'log')
%xlim([0, 1000])
ylim([0 0.015])

subplot(4,1,2)
bar(1:1:num_bcds, bcdabund(rankt2,2))
xlabel('barcode')
ylabel('abundance')
title('t=48h descending order')
%set(gca,'FontSize',20,'LineWidth',1.5)
set(gca,'FontSize',20,'LineWidth',1.5, 'YScale', 'log')
%xlim([0, 1000])
ylim([0 0.015])

subplot(4,1,3)
bar(1:1:num_bcds, bcdabund(rankt3,3))
xlabel('barcode')
ylabel('abundance')
title('t=96h descending order')
%set(gca,'FontSize',20,'LineWidth',1.5)
set(gca,'FontSize',20,'LineWidth',1.5,'YScale', 'log')
%xlim([0, 1000])
ylim([0 0.015])

subplot(4,1,4)
bar(1:1:num_bcds, bcdabund(rankt4,4))
xlabel('barcode')
ylabel('abundance')
title('t=144h descending order')
%set(gca,'FontSize',20,'LineWidth',1.5)
set(gca,'FontSize',20,'LineWidth',1.5,'YScale', 'log')
%xlim([0, 1000])
ylim([0 0.015])
edges = linspace(0,0.01, 500);
pmean = mean(bcdabund)
pstdev = std(bcdabund)

%%  Barcode Frequency histogram separated by timepoint
figure;
subplot(4,1,1)
histogram(bcdabund(:,1),edges)
xlabel('Abundance')
ylabel('Frequency')
title('t=0h')
%set(gca,'FontSize',20,'LineWidth',1.5,'XScale', 'log', 'YScale', 'log')
set(gca,'FontSize',20,'LineWidth',1.5)
xlim([0, 0.001])
ylim([0, 5000])

subplot(4,1,2)
histogram(bcdabund(:,2),edges)
xlabel('Abundance')
ylabel('Frequency')
title('t=48h')
%set(gca,'FontSize',20,'LineWidth',1.5,'XScale', 'log', 'YScale', 'log')
set(gca,'FontSize',20,'LineWidth',1.5)
xlim([0, 0.001])
ylim([0, 5000])

subplot(4,1,3)
histogram(bcdabund(:,3), edges)
xlabel('Abundance')
ylabel('Frequency')
title('t=96h')
%set(gca,'FontSize',20,'LineWidth',1.5,'XScale', 'log', 'YScale', 'log')
set(gca,'FontSize',20,'LineWidth',1.5)
xlim([0, 0.001])
ylim([0, 5000])

subplot(4,1,4)
histogram(bcdabund(:,4), edges)
xlabel('Abundance')
ylabel('Frequency')
title('t=144h')
%set(gca,'FontSize',20,'LineWidth',1.5,'XScale', 'log', 'YScale', 'log')
set(gca,'FontSize',20,'LineWidth',1.5)
xlim([0, 0.001])
ylim([0, 5000])
%% Barcode frequency histogram on one plot
figure;
histogram(bcdabund(:,1),edges, 'FaceAlpha', 0.5)
hold on
%set(gca,'FontSize',20,'LineWidth',1.5,'XScale', 'log', 'YScale', 'log')
%set(gca,'FontSize',20,'LineWidth',1.5)
xlim([0, 0.001])
%ylim([0, 5000])
histogram(bcdabund(:,2),edges, 'FaceAlpha', 0.5)
histogram(bcdabund(:,3), edges, 'FaceAlpha', 0.5)
histogram(bcdabund(:,4), edges, 'FaceAlpha', 0.5)
xlabel('Abundance')
ylabel('Frequency')
title('All time points overlayed')
legend('t=0h', 't=48h', 't=96h', 't=144h')
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5, 'YScale', 'log')

%% Calculate Shannon entropy
% H = sum(pi log pi)
p1 = bcdabund(:,1);
H1 = -sum(p1.*log2(p1));
p2 = bcdabund(:,2);
H2 = -sum(p2.*log2(p2));
p3 = bcdabund(:,3);
H3 = -sum(p3.*log2(p3));
p4 = bcdabund(:,4);
H4 = -sum(p4.*log2(p4));

Hvec = vertcat(H1, H2, H3, H4)

figure;
plot([0,48,96,144], Hvec, 'bo', 'LineWidth', 10)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('time (hours)')
ylabel('Shannon Entropy')
title('Shannon Entropy H at each time point')

% Simpson Index
D1 = 1./(sum(p1.^2));
D2 = 1./(sum(p2.^2));
D3 = 1./(sum(p3.^2));
D4 = 1./(sum(p4.^2));
Dvec = vertcat(D1, D2, D3, D4);

figure;
plot([0,48,96,144], Dvec, 'ro', 'LineWidth', 10)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('time (hours)')
ylabel('Simpson Index')
title('Simpson Index D at each time point')
%ylim([50, 900])

%% Diversity index calculation 
qvec = [];
%qvec = [1e-2, 1e-1, 1, 10, 100, 1000];
qvec = logspace(-2, 5, 20)
% qd = (sum(pi^q))^(1/(1-q))
Qd1 = [];
Qd2 = [];
Qd3 = [];
Qd4= [];
for i = 1:length(qvec)
Qd1(i) = (sum(p1.^qvec(i))).^(1./(1-qvec(i)));
Qd2(i)= (sum(p2.^qvec(i))).^(1./(1-qvec(i)));
Qd3(i)=(sum(p3.^qvec(i))).^(1./(1-qvec(i)));
Qd4(i) =(sum(p4.^qvec(i))).^(1./(1-qvec(i)));
end

figure;
plot(qvec, Qd1, 'o-', 'LineWidth', 3)
hold on
plot(qvec, Qd2, 'o-', 'LineWidth', 3)
plot(qvec, Qd3, 'o-', 'LineWidth', 3)
plot(qvec, Qd4, 'o-', 'LineWidth', 3)
legend('t=0h', 't=48h', 't=96h', 't=144h')
legend boxoff
xlabel('Diverisity parameter q')
ylabel('Diversity Index qD')
title('Diversity index qD changes over time')
set(gca,'FontSize',20,'LineWidth',1.5, 'XScale', 'log')
%xlim([qvec(1) qvec(20)])

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

%% Plot the indiviudal growth trajectories for all cells
figure;
for i = 1:1000%num_bcds
    tvec =  bcd(i).time;
    nvec =  bcd(i).ntot(1:4);
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
    nvec =  bcd(i).ntot(1:4);
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
    nvec =  bcd(i).ntot(1:4);
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

figure;
for i = 1:1000%num_bcds
    tvec =  bcd(i).time;
    nvec =  bcd(i).ntot(1:4);
    %tvec =  bcd(i).time;
    %nvec =  bcd(i).ntot;
    plot(tvec, nvec/nvec(1), '-', 'LineWidth', 1)
    set(gca,'FontSize',20,'LineWidth',1.5)
    hold on
    xlabel('time (hours)')
    ylabel('N(t)/N_{0}')
    title('N_{0} normalized growth trajectories')
    xlim([tvec(1) tvec(end)])
end


%% Plot the indiviudal growth trajectories for cells colored by CXCR4, CD18
figure;
for i = 1:1000%num_bcds
    tvec =  bcd(i).time;
    nvec =  bcd(i).ntot(1:4);
    %tvec =  bcd(i).time;
    %nvec =  bcd(i).ntot;
    if bcd(i).CD18pos == 1
    plot(tvec, log(nvec/nvec(1)), 'b-', 'LineWidth', 0.2)
    end
    if bcd(i).CXCR4pos == 1
        plot(tvec, log(nvec/nvec(1)), 'r-', 'LineWidth', 1)
    end
    if bcd(i).CXCR4pos==0 && bcd(i).CD18pos==0
        plot(tvec, log(nvec/nvec(1)), 'k-', 'LineWidth', 1)
    end
    
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
    nvec =  bcd(i).ntot(1:4);
    %tvec =  bcd(i).time;
    %nvec =  bcd(i).ntot;
    if bcd(i).CD18pos == 1
    plot(tvec, nvec/nvec(1), 'b-', 'LineWidth', 0.2)
    end
    if bcd(i).CXCR4pos == 1
        plot(tvec, nvec/nvec(1), 'r-', 'LineWidth', 1)
    end
    if bcd(i).CXCR4pos==0 && bcd(i).CD18pos==0
        plot(tvec, nvec/nvec(1), 'k-', 'LineWidth', 1)
    end
  
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
    nvec =  bcd(i).ntot(1:4);
    %tvec =  bcd(i).time;
    %nvec =  bcd(i).ntot;
    if bcd(i).CD18pos == 1
    plot(tvec, nvec/nvec(1), 'b-', 'LineWidth', 0.2)
    end
    if bcd(i).CXCR4pos == 1
        plot(tvec, nvec/nvec(1), 'r-', 'LineWidth', 1)
    end
    if bcd(i).CXCR4pos==0 && bcd(i).CD18pos==0
        plot(tvec, nvec/nvec(1), 'k-', 'LineWidth', 1)
    end
  
    set(gca,'FontSize',20,'LineWidth',1.5)
    hold on
    xlabel('time (hours)')
    ylabel('N(t)')
    title('N_{0} normalized growth trajectories')
    xlim([tvec(1) tvec(end)])
end

%% Fit the trajectories from the 3 observed time points for each lineage



for i = 1:num_bcds
      tvec = bcd(i).time; % shift time points for fitting set t=48h to t=0 for now
      nvec = bcd(i).ntot(2:4);
      num_params = 1;
      % Fit on all but t=48h, which we will call t=0
      N0 = bcd(i).ntot(1); % set N0 as the first observed tp at t=48h
      nbar = mean(nvec); % average value of volumes over all measured times
      
      % First fit each curve to single exponential over entire time course
      LB = -Inf ;  % Lower Bounds
      UB = Inf; % Upper Bounds
      params0 = log(nvec(3)/nvec(2))./(48); % best guess (slope of the last two tps)
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
  
end
%% Plot some example results of fit versus data
figure;
for i = 1:10
    plot(bcd(i).time, bcd(i).ntot(1:4), '*')
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
greslist = [];
gsenslist = [];
idCXCR4pos = [];
for i = 1:num_bcds
    if  bcd(i).badfit ==0
    glist = vertcat(glist, bcd(i).g);
    
    if bcd(i).CXCR4pos==1
        greslist = vertcat(greslist, bcd(i).g);
        idCXCR4pos = vertcat(idCXCR4pos, 1);
    end
    if bcd(i).CXCR4pos == 0
        gsenslist = vertcat(gsenslist, bcd(i).g);
        idCXCR4pos = vertcat(idCXCR4pos, 0);
    end
    end    
end
[gdistrib, rankg] = sort(glist, 'descend');
gmean = mean(glist);
gmed = median(glist);
gstdev = std(glist);

[gresdistrib, rankgres] = sort(greslist, 'descend');
gresmean = mean(greslist);
gresmed = median(greslist);
gresstdev = std(greslist);

[gsensdistrib, rankgsens] = sort(gsenslist, 'descend');
gsensmean = mean(gsenslist);
gsensmed = median(gsenslist);
gsensstdev = std(gsenslist);

%% Plot growth rate distribution
figure;
bar(1:1:length(glist), glist)
xlabel('barcode')
ylabel('growth rate')
title('growth rate distribution')
set(gca,'FontSize',20,'LineWidth',1.5)
%xlim([0, 500])
%ylim([0 0.097])
%% Ordered growth rate distribution
figure;
bar(1:1:length(glist), gdistrib)
hold on
xlabel('barcode')
ylabel('growth rate')
title('growth rate distribution')
set(gca,'FontSize',20,'LineWidth',1.5)
%xlim([0, 500])
%ylim([0 0.097])

% Colored by CXCR4 pos/negative lineages
ordered_IDS = idCXCR4pos(rankg);

index = 1:1:length(gdistrib);
gdistribpos = gdistrib(ordered_IDS==1);
barcodespos = index(ordered_IDS==1);
gdistribneg = gdistrib(ordered_IDS==0);
barcodesneg = index(ordered_IDS ==0);
figure;
bar(barcodespos, gdistribpos, 'r')
hold on
bar(barcodesneg, gdistribneg, 'b', 'FaceAlpha', 0.4)
xlabel ('barcodes')
ylabel('Fitted growth rate')
legend('CXCR4+ lineages', 'CXCR4- lineages')
legend boxoff
title('Growth rate distribution stratified by CXCR4')
set(gca,'FontSize',20,'LineWidth',1.5)
%% Stratified dot plot
figure;
for i = 1:length(gdistrib)
    hold on
    if ordered_IDS(i)==0
    plot(ordered_IDS(i), gdistrib(i), 'b.')
    end
    if ordered_IDS(i)==1
        plot(ordered_IDS(i), gdistrib(i), 'r.')
    end
end
xlim([-0.5, 1.5])
ylabel('growth rate')
set(gca,'FontSize',20,'LineWidth',1.5,'XTickLabel',{'CXCR4-','CXCR4+'})
%%
group = [ 0*ones(size(gdistribneg)); ones(size(gdistribpos))];
figure;
boxplot([gdistribneg; gdistribpos],group)
set(gca,'FontSize',20,'LineWidth',1.5,'XTickLabel',{'CXCR4-','CXCR4+'})
ylabel('growth rate')


%%
figure;
boxplot([gdistribneg, gdistribpos], [zeros(size(gistribneg)), ones(size(gdistribpos))])
%%
figure;
histogram(greslist, 50,'Normalization', 'probability', 'FaceColor', 'r')
hold on
histogram(glist, 50,'Normalization', 'probability')
legend('CXCR4','population')
legend boxoff
xlabel('growth rate')
ylabel('frequency')
title('growth rate distribution')
set(gca,'FontSize',20,'LineWidth',1.5)

figure;
histogram(greslist,50, 'Normalization', 'probability', 'FaceColor', 'r')
hold on
histogram(gsenslist,50, 'Normalization', 'probability', 'FaceColor', 'b', 'FaceAlpha', 0.4)
legend('CXCR4+', 'CXCR4-')
legend boxoff
xlabel('growth rate')
ylabel('frequency')
title('Growth rate distributions CXCR4+ & CXCR4-')
set(gca,'FontSize',20,'LineWidth',1.5)
%% Look at the number of barcodes vs. time barcodes observed

%Make a frequency table
%first column = times barcode observed
%second column = number of barcodes
%third column = percent of total barcodes
ftable = [];

ftable0= tabulate(Nall(:,1));
ftable0(:,3)=ftable0(:,3)./100;
ftable48= tabulate(Nall(:,2));
ftable48(:,3)=ftable48(:,3)./100;
ftable96 = tabulate(Nall(:,3));
ftable96(:,3) = ftable96(:,3)./100;
ftable144 = tabulate(Nall(:,4));
ftable144(:,3) = ftable144(:,3)./100;
%%

figure;
plot(ftable0(:,1), ftable0(:,2), '.')
hold on
plot(ftable48(:,1), ftable48(:,2), '.')
plot(ftable96(:,1), ftable96(:,2), '.')
plot(ftable144(:,1), ftable144(:,2), '.')
legend('t=0h', 't=48h', 't=96h', 't=144h')
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5, 'Xscale', 'log', 'Yscale', 'log')
xlabel ('Times barcode is observed')
ylabel('Number of barcodes')
%%
figure;
plot(ftable0(:,1), ftable0(:,3), '.')
hold on
plot(ftable48(:,1), ftable48(:,3), '.')
plot(ftable96(:,1), ftable96(:,3), '.')
plot(ftable144(:,1), ftable144(:,3), '.')
legend('t=0h', 't=48h', 't=96h', 't=144h')
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5, 'Xscale', 'log', 'Yscale', 'log')
xlabel ('Times barcode is observed')
ylabel('Relative abundance of barcode')
