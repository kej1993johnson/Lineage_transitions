% Barcode sampling 
% This script loads in the barcode sampling data which contains raw numbers
% of barcodes found in samples at three different time points after initial
% barcoding:
%OGTP1 = 48h, 3% population which was 19.8x10^6 cells
%OGTP2 = 96h, 3% population which was 75.9x10^6 cells 
%OGTP3 = 144h, 3% population was at 180.8x10^6 cells

% We will reset these times to t=0, 48, and 96 h for modeling purposes. 

% For the CD18 and CXCR4 populations the number of cells sampled is:
% CD18 cells = 431,200 % sensitive
% CXCR4 cells = 71,867 % resistant

close all; clear all; clc
%% Load in OG and TPO barcode frequency distribution
[N, T]= xlsread('../data/Barcode_sampling_pre_treat.xlsx');
% Want to make a structure that tracks the number of barcodes, the barcode
% the rank based on the first time point, second time point, and last time
% point and the scaled population number


nreads = nansum(N); % sum each column, to get total number of reads before 
% removing.

% This contains columns of data that consist of 
% t=48h, t=96h, t=144h, CD18+, CXCR4+
%% Remove barcodes not observed in all the first time point
% Want to make a structure that tracks the number of barcodes, the barcode
% the rank based on the first time point, second time point, and last time
% point and the scaled population number

% Want to set all NaNs to 0 in the N matrix
idnotobs = isnan(N(:,1)); % We want only barcodes observed in all three time points 
Ncons=[];
cons_barcodes = {};
for i = 1:length(N)
    if ~any(idnotobs(i,:))
        Ncons = vertcat(Ncons, N(i,:));
        cons_barcodes = vertcat(cons_barcodes, char(T(i+1,1)));
    end
end
% Result is we go from 14004 barcodes to 8107 consensus barcodes
nreadscons = nansum(Ncons);
N(isnan(N))=0;
Ncons(isnan(Ncons))=0;
lostreads = nreads-nreadscons;
proplost = lostreads./nreads % proportion lost increases 

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
numsampscons = numsamps.*(1-proplost);
numtotscons = numtots.*(1-proplost); % look only at the total population size that we have consensus barcodes for
% Note that for the numtots for the CD18 and CXCR4 populations this may not
% be right, but since we're not using them for growth analysis it's ok.
pct_samps = numsamps(1:3)./numtots(1:3); 
 % these should all be around 3 but just to check

% To estimate the number of cells in each barcode in the sample and the
% population you need
bcdabund = N./nreads; % proportion of barcodes at each time point
bcdabundcons = Ncons./nreadscons;
Nsamps = bcdabundcons.*numsampscons; % number of cells with that barcode in the sample
Ntot = bcdabundcons.*numtotscons; % number of cells with that barcode at each time point

%% Check that the concensus barcode abundances add to 1
probtots = sum(bcdabundcons(:,1:3))
% ensure that these represent real pdfs
%%
bcd = struct( 'time', zeros(1,1),'barcode',strings(1), 'nsamp',...
    zeros(1,1),'ntot', zeros(1,1),'abund', zeros(1,1), ...
    'CD18pos', zeros(1,1), 'CXCR4pos', zeros(1,1),'rank1',...
    zeros(1,1), 'rank2', zeros(1,1), 'rank3', zeros(1,1));
% Make the rankt1, rankt2, and rankt3 list of indices
% Use these to keep track of how the barcode moves in ranking
[B1, rankt1]=sort(Ncons(:,1), 'descend');
[B2, rankt2] = sort(Ncons(:,2), 'descend');
[B3, rankt3] = sort(Ncons(:,3), 'descend');
%% Loop through and add each value to the barcode structure
num_bcds = length(Ncons);
for i = 1:num_bcds
    bcd(i).time = [0, 48, 96]; % RESET THE TIME SO THAT 48h from thaw is t=0h
    bcd(i).barcode = char(T(i+1,1));
    bcd(i).nsamp = Nsamps(i,:);
    bcd(i).ntot = Ntot(i,:);
    bcd(i).abund = bcdabundcons(i,:);
    bcd(i).rank1 = rankt1(i);
    bcd(i).rank2 =rankt2(i);
    bcd(i).rank3 = rankt3(i);
end
%% Look at the barcode distribution histogram, here ordered by t=0 barcodes
figure;
subplot(3,1,1)
bar(1:1:num_bcds, (bcdabundcons(:,1)))
xlabel('barcode')
ylabel('abundance')
title('t=0h barcode abundance distribution')
ylim([0 0.1])
xlim([0, num_bcds])
set(gca,'FontSize',20,'LineWidth',1.5, 'YScale', 'log')

subplot(3,1,2)
bar(1:1:num_bcds, (bcdabundcons(:,2)))
xlabel('barcode')
ylabel('abundance')
title('t=48h ordered by t=0h')
ylim([0 0.1])
xlim([0, num_bcds])
set(gca,'FontSize',20,'LineWidth',1.5, 'YScale', 'log')

subplot(3,1,3)
bar(1:1:num_bcds, (bcdabundcons(:,3)))
xlabel('barcode')
ylabel('abundance')
title('t=96h ordered by t=0h')
ylim([0 0.1])
xlim([0, num_bcds])
set(gca,'FontSize',20,'LineWidth',1.5, 'YScale', 'log')
%% Order by abundance
figure;
subplot(3,1,1)
bar(1:1:num_bcds, (bcdabundcons(:,1)))
xlabel('barcode')
ylabel('abundance')
title('t=0h barcode abundance distribution')
ylim([0 0.1])
xlim([0, num_bcds])
set(gca,'FontSize',20,'LineWidth',1.5, 'YScale', 'log')

subplot(3,1,2)
bar(1:1:num_bcds, (bcdabundcons(rankt2,2)))
xlabel('barcode')
ylabel('abundance')
title('t=48h ordered by abundance')
ylim([0 0.1])
xlim([0, num_bcds])
set(gca,'FontSize',20,'LineWidth',1.5, 'YScale', 'log')

subplot(3,1,3)
bar(1:1:num_bcds, (bcdabundcons(rankt3,3)))
xlabel('barcode')
ylabel('abundance')
title('t=96h ordered by abundance')
ylim([0 0.1])
xlim([0, num_bcds])
set(gca,'FontSize',20,'LineWidth',1.5, 'YScale', 'log')
%% Overlaid

figure;
bar(1:1:num_bcds, (bcdabundcons(:,1)))
hold on
bar(1:1:num_bcds, (bcdabundcons(rankt2,2)))
bar(1:1:num_bcds, (bcdabundcons(rankt3,3)))
xlabel('barcode')
ylabel('abundance')
title('Barcode Abundance at each time point')
legend('t=0h', 't=48h', 't=96h')
legend boxoff
ylim([0 0.1])
xlim([0, num_bcds])
set(gca,'FontSize',20,'LineWidth',1.5, 'YScale', 'log')

%%


% make this into a 3 column heat map with barcode on the y-axis and time
% point on the x axis, colored by bcdabundance.
xlabs = {'t=0h', 't=48h', 't=96h'};
figure;
colormap('hot')
imagesc(log10(bcdabundcons(:, 1:3)))
h = colorbar;
h.Label.String = 'Log10 barcode abundance'
ylabel(h, 'Log10 barcode abundance')
set(gca,'FontSize',16,'LineWidth',1.5)
ylabel('barcode')
title('Heatmap of barcode abundance')
set(gca,'XTick',[1, 2 ,3], 'XTickLabel', xlabs);% something like this