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
%% Remove barcodes not observed in all 3 time points
% Want to make a structure that tracks the number of barcodes, the barcode
% the rank based on the first time point, second time point, and last time
% point and the scaled population number

% Want to set all NaNs to 0 in the N matrix
idnotobs = isnan(N(:,1:3)); % We want only barcodes observed in all three time points 
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
%% Export the population size for each lineage as a csv for reading in EvoFreq
smallNtot = Ntot;
%smallNtot = Ntot/1e4;
parentsw = 1:1:length(Ntot);
parents = parentsw';
clones = parentsw';
t1 = round(smallNtot(:,1),0);
t2 = round(smallNtot(:,2),0);
t3 = round(smallNtot(:,3),0);

Tevowide = table(parents, clones, t1, t2, t3)

writetable(Tevowide, 'Tdata.csv')

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

%% Make a list of the barcodes ordered by their magnitude in the third time point
cons_barcodes1 = cons_barcodes(rankt1);
cons_barcodes2 = cons_barcodes(rankt2);
cons_barcodes3 = cons_barcodes(rankt3);
Tbcds = table(cons_barcodes1, cons_barcodes2, cons_barcodes3);

writetable(Tbcds, 'ordered_barcodes.csv')
%% Loop through and add each value to the barcode structure
num_bcds = length(Ncons);
thres_ratio = 1; % abundance to call a lineage CD18 positive
CD18ratio = [];
ct_CXCR4pos = 0;
for i = 1:num_bcds
    bcd(i).time = [0, 48, 96]; % RESET THE TIME SO THAT 48h from thaw is t=0h
    bcd(i).barcode = char(T(i+1,1));
    bcd(i).nsamp = Nsamps(i,:);
    bcd(i).ntot = Ntot(i,:);
    bcd(i).abund = bcdabundcons(i,:);
    bcd(i).rank1 = rankt1(i);
    bcd(i).rank2 =rankt2(i);
    bcd(i).rank3 = rankt3(i);
    bcd(i).CD18ratio = bcd(i).abund(1,4)/bcd(i).abund(1,5); % CD18 abund/ CXCR4 abund
    CD18ratio = vertcat(CD18ratio, bcd(i).CD18ratio);
    if bcd(i).CD18ratio >1
        bcd(i).CD18pos = 1;
        bcd(i).CXCR4pos = 0;
    end
    if bcd(i).CD18ratio <1
        bcd(i).CD18pos = 0;
        bcd(i).CXCR4pos = 1;
        ct_CXCR4pos= ct_CXCR4pos +1;
    end
    if isnan(bcd(i).CD18ratio)
        bcd(i).CD18pos = 0;
        bcd(i).CXCR4pos = 0;
    end
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
%%  Barcode Frequency histogram separated by timepoint
edges = linspace(0,0.01, 500);
figure;
subplot(3,1,1)
histogram(bcdabundcons(:,1),edges)
xlabel('Abundance')
ylabel('Frequency')
title('t=0h')
%set(gca,'FontSize',20,'LineWidth',1.5,'XScale', 'log', 'YScale', 'log')
set(gca,'FontSize',20,'LineWidth',1.5)
xlim([0, 0.001])
ylim([0, 5000])

subplot(3,1,2)
histogram(bcdabundcons(:,2),edges)
xlabel('Abundance')
ylabel('Frequency')
title('t=48h')
%set(gca,'FontSize',20,'LineWidth',1.5,'XScale', 'log', 'YScale', 'log')
set(gca,'FontSize',20,'LineWidth',1.5)
xlim([0, 0.001])
ylim([0, 5000])

subplot(3,1,3)
histogram(bcdabundcons(:,3), edges)
xlabel('Abundance')
ylabel('Frequency')
title('t=96h')
%set(gca,'FontSize',20,'LineWidth',1.5,'XScale', 'log', 'YScale', 'log')
set(gca,'FontSize',20,'LineWidth',1.5)
xlim([0, 0.001])
ylim([0, 5000])

% Barcode frequency histogram on one plot
figure;
histogram(bcdabundcons(:,1),edges, 'FaceAlpha', 0.5)
hold on
%set(gca,'FontSize',20,'LineWidth',1.5,'XScale', 'log', 'YScale', 'log')
%set(gca,'FontSize',20,'LineWidth',1.5)
xlim([0, 0.001])
%ylim([0, 5000])
histogram(bcdabundcons(:,2),edges, 'FaceAlpha', 0.5)
histogram(bcdabund(:,3), edges, 'FaceAlpha', 0.5)
xlabel('Abundance')
ylabel('Frequency')
title('All time points overlayed')
legend('t=0h', 't=48h', 't=96h')
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5)% 'YScale', 'log')
%% Calculate Shannon entropy
% H = sum(pi log pi)
p1 = bcdabundcons(:,1);
H1 = -sum(p1.*log2(p1));
p2 = bcdabundcons(:,2);
H2 = -sum(p2.*log2(p2));
p3 = bcdabundcons(:,3);
H3 = -sum(p3.*log2(p3));

Hvec = vertcat(H1, H2, H3);

figure;
plot([0,48,96], Hvec, 'bo', 'LineWidth', 10)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('time (hours)')
ylabel('Shannon Entropy')
title('Shannon Entropy H at each time point')

% Simpson Index
D1 = 1./(sum(p1.^2));
D2 = 1./(sum(p2.^2));
D3 = 1./(sum(p3.^2));
Dvec = vertcat(D1, D2, D3);

figure;
plot([0,48,96], Dvec, 'ro', 'LineWidth', 10)
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
for i = 1:length(qvec)
Qd1(i) = (sum(p1.^qvec(i))).^(1./(1-qvec(i)));
Qd2(i)= (sum(p2.^qvec(i))).^(1./(1-qvec(i)));
Qd3(i)=(sum(p3.^qvec(i))).^(1./(1-qvec(i)));
end

figure;
plot(qvec, Qd1, 'o-', 'LineWidth', 3)
hold on
plot(qvec, Qd2, 'o-', 'LineWidth', 3)
plot(qvec, Qd3, 'o-', 'LineWidth', 3)
legend('t=0h', 't=48h', 't=96h')
legend boxoff
xlabel('Diverisity parameter q')
ylabel('Diversity Index qD')
title('Diversity index qD changes over time')
set(gca,'FontSize',20,'LineWidth',1.5, 'XScale', 'log')
%xlim([qvec(1) qvec(20)])
%% Plot lineage growth trajectories colored by CXCR4 +/-
figure;
for i = 1:1000%num_bcds
    tvec =  bcd(i).time;
    nvec =  bcd(i).ntot(1:3);
    %tvec =  bcd(i).time;
    %nvec =  bcd(i).ntot;
    if bcd(i).CD18pos == 1
    plot(tvec, nvec/nvec(1), 'b-', 'LineWidth', 0.1)
    end
    if bcd(i).CXCR4pos == 1
        plot(tvec, nvec/nvec(1), 'r-', 'LineWidth', 2)
    end
    if bcd(i).CXCR4pos==0 && bcd(i).CD18pos==0
        plot(tvec, nvec/nvec(1), 'k-', 'LineWidth', 1)
    end
  
    set(gca,'FontSize',20,'LineWidth',1.5)
    hold on
    xlabel('time (hours)')
    ylabel('N(t)')
    title('N_{0} normalized lineage growth rates')
    xlim([tvec(1) tvec(end)])
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



%% Barcode abundance distribution reordered in each time point
figure;
subplot(3,1,1)
bar(1:1:num_bcds, bcdabundcons(:,1))
xlabel('barcode')
ylabel('abundance')
title('t=0h barcode abundance distribution')
set(gca,'FontSize',20,'LineWidth',1.5, 'YScale', 'log')
%xlim([0, 1000])
ylim([0 0.1])

subplot(3,1,2)
bar(1:1:num_bcds, bcdabundcons(rankt2,2))
xlabel('barcode')
ylabel('abundance')
title('t=48h descending order')
set(gca,'FontSize',20,'LineWidth',1.5, 'YScale', 'log')
%xlim([0, 1000])
ylim([0 0.1])

subplot(3,1,3)
bar(1:1:num_bcds, bcdabundcons(rankt3,3))
xlabel('barcode')
ylabel('abundance')
title('t=96h descending order')
set(gca,'FontSize',20,'LineWidth',1.5,'YScale', 'log')
%xlim([0, 1000])
ylim([0 0.1])
%% Plot as separate figures
figure;
bar(1:1:num_bcds, bcdabundcons(:,1))
xlabel('barcode')
ylabel('abundance')
title('t=0h')
set(gca,'FontSize',20,'LineWidth',1.5,'YScale', 'log')
%xlim([0, 500])
ylim([0 0.097])

figure;
bar(1:1:num_bcds, bcdabundcons(rankt2,2))
xlabel('barcode')
ylabel('abundance')
title('t=48h ')
set(gca,'FontSize',20,'LineWidth',1.5, 'YScale', 'log')
%xlim([0, 500])
ylim([0 0.097])

figure;
bar(1:1:num_bcds, bcdabundcons(rankt3,3))
xlabel('barcode')
ylabel('abundance')
title('t=96h')
set(gca,'FontSize',20,'LineWidth',1.5, 'YScale', 'log')
%xlim([0, 500])
ylim([0 0.097])

%% Plot the indiviudal growth trajectories
figure;
for i = 1:1000%num_bcds
    tvec =  bcd(i).time;
    nvec =  bcd(i).ntot(1:3);
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
%%
figure;
for i = 1:1000%num_bcds
    tvec =  bcd(i).time;
    nvec =  bcd(i).ntot(1:3);
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
%%
figure;
for i = 1:1000%num_bcds
    tvec =  bcd(i).time;
    nvec =  bcd(i).ntot(1:3);
    %tvec =  bcd(i).time;
    %nvec =  bcd(i).ntot;
    plot(tvec, nvec/nvec(1), '-', 'LineWidth', 1)
    set(gca,'FontSize',20,'LineWidth',1.5)
    hold on
    xlabel('time (hours)')
    ylabel('N(t)')
    title('N_{0} Normalized lineage growth trajectories')
    xlim([tvec(1) tvec(end)])
end
%% Fit the trajectories from the 3 observed time points for each lineage



for i = 1:num_bcds
      tvec = bcd(i).time; % shift time points for fitting set t=48h to t=0 for now
      nvec = bcd(i).ntot(2:3);
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
end
%% Plot some example results of fit versus data
figure;
for i = 1:10
    plot(bcd(i).time, bcd(i).ntot(1:3), '*')
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
legend('CXCR4+ lineages', 'CD18+ lineages')
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
set(gca,'FontSize',20,'LineWidth',1.5,'XTickLabel',{'CD18+','CXCR4+'})
%%
group = [ 0*ones(size(gdistribneg)); ones(size(gdistribpos))];
figure;
boxplot([gdistribneg; gdistribpos],group)
set(gca,'FontSize',20,'LineWidth',1.5,'XTickLabel',{'CD18+','CXCR4+'})
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
edgeshist = linspace(-0.02,0.1, 100);
figure;
histogram(greslist,edgeshist, 'Normalization', 'probability', 'FaceColor', 'r')
hold on
histogram(gsenslist, edgeshist, 'Normalization', 'probability', 'FaceColor', 'b', 'FaceAlpha', 0.4)
legend('CXCR4+', 'CD18+')
legend boxoff
xlabel('growth rate')
ylabel('frequency')
title('Growth rate distributions CXCR4+ & CD18+')
set(gca,'FontSize',20,'LineWidth',1.5)


gresmean = mean(greslist)
gresmed = median(greslist)
gresstdev = std(greslist)
gsensmean = mean(gsenslist)
gsensmed = median(gsenslist)
gsensstdev = std(gsenslist)


%% Plot growth rate distribution
glist = [];
for i = 1:num_bcds
    if bcd(i).badfit ==0
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

%% Look at the number of barcodes vs. time barcodes observed

%Make a frequency table
%first column = times barcode observed
%second column = number of barcodes
%third column = percent of total barcodes
ftable = [];

ftable0= tabulate(Ncons(:,1));
ftable0(:,3)=ftable0(:,3)./100;
ftable48= tabulate(Ncons(:,2));
ftable48(:,3)=ftable48(:,3)./100;
ftable96 = tabulate(Ncons(:,3));
ftable96(:,3) = ftable96(:,3)./100;
%%

figure;
plot(ftable0(:,1), ftable0(:,2), '.')
hold on
plot(ftable48(:,1), ftable48(:,2), '.')
plot(ftable96(:,1), ftable96(:,2), '.')
legend('t=0h', 't=48h', 't=96h')
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5, 'Xscale', 'log', 'Yscale', 'log')
xlabel ('Times barcode is observed')
ylabel('Number of barcodes')

figure;
plot(ftable0(:,1), ftable0(:,3), '.')
hold on
plot(ftable48(:,1), ftable48(:,3), '.')
plot(ftable96(:,1), ftable96(:,3), '.')
legend('t=0h', 't=48h', 't=96h', 't=144h')
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5, 'Xscale', 'log', 'Yscale', 'log')
xlabel ('Times barcode is observed')
ylabel('Relative abundance of barcode')

%% Make stacked bar graph of proportion in CXCR4+ CD18+ cells over time.
CXCR4abund_t= [0 0 0; 0 0 0 ];
xlabs = {'t=0h', 't=48h', 't=96h'}
for i = 1:num_bcds
    if bcd(i).CXCR4pos==1
    CXCR4abund_t(2,:) = CXCR4abund_t(2,:) + bcd(i).abund(1:3);
    end
end
CXCR4abund_t(1,:) = 1-CXCR4abund_t(2,:);
figure;
bar(CXCR4abund_t', 'stacked')
set(gca, 'FontSize', 14', 'XTickLabel', xlabs)
legend('CD18+', 'CXCR4+')
%%
figure;
plot([ 0 48 96], CXCR4abund_t(2,:), 'r*-','LineWidth', 3)
hold on
plot([ 0 48 96], 1-CXCR4abund_t(2,:),'b*-', 'LineWidth', 3)
set(gca, 'FontSize', 14)
xlabel ('Time (hours)')
ylabel('Abundance')
legend('CXCR4+', 'CD18+')
legend boxoff
ylim([ 0 1.1])
xlim([-10 106])
