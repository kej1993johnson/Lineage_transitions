% Model of growth rate distributions and the effect on barcoding frequency
close all; clear all; clc
% Find apparent growth rate
% if at t=0d, sort out 5000 cells, for an initial diversity of 5000, then
% apparent growth rate from 
% N(t=22d) = N(t=0)e^(g*t)
% Nt/N0= e^(g*t)
% ln both sides
% ln(Nt)-ln(N0) = g*t
% g = (ln(Nt)-ln(N0))/(t)
Nts=3.382e12; % N at time of scRNA seq (completely arbitrary guess)
N0=5000; % number of barcoded cells sorted initially
tgr=22; % time of outgrowth before scRNAseq
g = (log(Nts)-log(N0))/tgr; % in units of per cell per day

% Create a normal distribution around this growth rate with an arbitrary
% standard deviation
stdev= 0.05;
gvals = normrnd(g,stdev,100);

% make histogram to visualize
histogram(gvals)
xlabel('growth rate')
xlim([0.5 1.3])
ylabel('frequency')
title('Normal growth rate distribution, g = 0.92 cells/day, \sigma = 0.05')

cell_growth_rates = normrnd(g,stdev,N0,1);
%%
N_cells = ones(N0,1);
for j = 1:length(cell_growth_rates)
    N_cells(j,2)= exp(cell_growth_rates(j).*tgr);
end

sum_Ncells= sum(N_cells);
N0pop=sum_Ncells(1);
Nfpop=sum_Ncells(2);
apparent_g = (log(Nfpop)-log(N0pop))/tgr;
for j = 1:N0
sim_pct_abund(j) = 100*(N_cells(j,2)./Nfpop);
end
sim_bcd_abund= sort(sim_pct_abund, 'descend');
%
figure;
hold off
bar(1:1:50, sim_bcd_abund(1:50))
xlim([0 50])
ylabel('Percent Abundance')
xlabel('Barcode')
title('Simulated barcode abundance for normally distributed growth rates')


%% Now try back calculating from barcode frequencies 
% If we assume a doubling time of 18 hours for the cell population
% First find the total number of cells expected at 22 days
doubling_t = 18/24;
doubles = tgr/doubling_t;
N_tgr= 2^doubles*N0;
g_Hc3s = (log(N_tgr)-log(N0))/tgr;
x = 0:1:5000;
% generate a simulated plot of number of cells vs individul barcode
sim_barcode_freq = 6500*pdf('Exponential',x, 35);
barcode_abund = round(sim_barcode_freq,0);

figure;
bar(x, barcode_abund)
xlim([0,50])
xlabel('barcode')
ylabel('Number of cells with barcode')
title('Simulated barcode frequency')

% for each x, find it's percent abundance
Nsc= sum(barcode_abund); % in reality, this would be from 10x data total number of cells sampled
perc_abund =100*(barcode_abund/Nsc);

figure;
bar(x, perc_abund)
xlim([0,50])
xlabel('barcode')
ylabel('Percent of total cells')
title('Simulated barcode percent')

% Find the expected cell numbers based on the barcoded frequencies
cell_num_by_barcode = N_tgr.*(perc_abund/100);
g_by_barcode= log(cell_num_by_barcode)/tgr;

figure;
histogram(g_by_barcode, 100)
xlabel('Growth rate')
ylabel('Frequency')
title('Simulated growth rate distribution')