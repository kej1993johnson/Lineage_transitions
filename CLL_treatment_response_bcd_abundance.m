% Analyze barcode abundance data taken at time of treatment response
% (trough)
% Cathy set up the experiment as follows: She plated 4 replicate of 10e6
% cells in 4 separate flasks, treated, and barcode sampled at Day 12
% post-treatment (need to confirm what the treatmen was). The number of
% reads for each barcode are in each of the columns

close all; clear all; clc;

%% We will use the old barcode sampling data to identify which lineages are
% CXCR4+ and which are CD8+

[N, T]= xlsread('../data/Barcode_sampling_pre_treat.xlsx');
N(isnan(N))=0;
%[T12]= importdata('../data/T50_barcode_analysis.txt')
nreads = nansum(N); 
bcdabund = N./nreads; %
CD18ratio = bcdabund(:,4)./bcdabund(:,5);
for i = 1:length(N)
bcd(i).barcode = char(T(i+1,1));
if CD18ratio(i)>1
bcd(i).CD18pos = 1;
bcd(i).CXCR4pos = 0;
end
if CD18ratio(i)<=1
    bcd(i).CD18pos = 0;
    bcd(i).CXCR4pos = 1;
end
end

CXCR4hibcds = [];
CD18hibcds = [];

for i = 1:length(bcd)
if bcd(i).CXCR4pos==1
    CXCR4hibcds = vertcat(CXCR4hibcds, bcd(i).barcode);
end
if bcd(i).CD18pos==1
    CD18hibcds = vertcat(CD18hibcds, bcd(i).barcode);
end

end

TCXCR4= table(CXCR4hibcds);
TCD18=table(CD18hibcds);

writetable(TCXCR4, '../out/barcodes_CXCR4_hi.csv')
writetable(TCD18, '../out/barcodes_CD18_hi.csv')
%% Check to see if the barcodes in the 12 DPT experiment are found in the list of CD18 and 
% CXCR4 positive barcodes
[N12, T12]= xlsread('../data/post_treat_bcd_abundance.xlsx');
[N0, T0] = xlsread('../data/TP0_lins.xlsx');
N12(isnan(N12))=0;
N0sum = N0(:,1)+N0(:,2);
[N0sumord, ordN0]=sort(N0sum, 'descend');
bcdsT0 = T0(2:end,1);
[N50, T50] = xlsread('../data/TP50_lins.xlsx');
N50sum = N50(:,1)+N50(:,2);
[N50sumord, ordN50]=sort(N50sum, 'descend');
bcdsT50=T50(2:end,1);

% Break this up into one structure with 5 separate matrices
bcd_treat(1).rep='10x Replicate';
bcd_treat(1).barcodes = T12(2:end,1);
bcd_treat(1).Nreads = N12(:,1);
bcd_treat(1).totreads= sum(N12(:,1));
bcd_treat(1).abund = bcd_treat(1).Nreads./bcd_treat(1).totreads;
bcd_treat(2).rep = 'Replicate 1';
bcd_treat(3).rep = 'Replicate 2';
bcd_treat(4).rep = 'Replicate 3';
bcd_treat(5).rep = 'Replicate 4';
bcd_treat(6).rep = 'TP0 (10x)';
bcd_treat(7).rep = 'TP1(10x)';

for i = 1:4
bcd_treat(i+1).barcodes = T12(2:nnz(N12(:,2*i+1))+1,2*i+1);
bcd_treat(i+1).Nreads = N12(1:nnz(N12(:,2*i+1)), 3);
bcd_treat(i+1).totreads = sum(bcd_treat(i+1).Nreads);
bcd_treat(i+1).abund = bcd_treat(i+1).Nreads./bcd_treat(i+1).totreads;
end

bcd_treat(6).barcodes = bcdsT0(ordN0,1);
bcd_treat(6).Nreads = N0sum(ordN0); % sum up the cells in high and low tol
bcd_treat(6).totreads= sum(N0sum); % tot number of cells
bcd_treat(6).abund = bcd_treat(6).Nreads./bcd_treat(6).totreads;
bcd_treat(7).barcodes = bcdsT50(ordN50,1);
bcd_treat(7).Nreads = N50sumord; % sum up the cells in high and low tol
bcd_treat(7).totreads= sum(N50sum); % tot number of cells
bcd_treat(7).abund = bcd_treat(7).Nreads./bcd_treat(7).totreads;
%% Test out a barcode abundance plot
figure;
bar(1:1:length(bcd_treat(6).abund), bcd_treat(6).abund)
xlabel('barcode')
ylabel('Abundance')
title('TP0(10x)')
set(gca,'FontSize',20,'LineWidth',1.5)
%% Go through each barcode and label by CXCR4 status
for j = 1:length(bcd_treat)
    % default to 0 in status
    CXCR4vec = zeros(length(bcd_treat(j).barcodes),1);
    CD18vec = zeros(length(bcd_treat(j).barcodes),1);
    % if barcode matches one of the list of CXCR4/CD18 barcodes, label
    for i = 1:length(bcd_treat(j).barcodes)
        if any(strcmp(TCXCR4.CXCR4hibcds,bcd_treat(j).barcodes(i)))
            CXCR4vec(i)=1;
        end
        if any(strcmp(TCD18.CD18hibcds,bcd_treat(j).barcodes(i)))
            CD18vec(i)=1;
        end
    end
    bcd_treat(j).CXCR4vec = CXCR4vec;
    bcd_treat(j).CD18vec = CD18vec;
end
%% Try plotting barcode abundance colored by CXCR4 and CD18 status
%figure;
for j = 1:length(bcd_treat)
    %subplot(2,3,j)
    figure;
  index = [];
  abundCXCR4hi = [];
  indCXCR4hi = [];
  abundCD18hi = [];
  indCD18hi = [];
  indneut = [];
  abundneut = [];
 index = 1:1:length(bcd_treat(j).abund);
abundCXCR4hi = bcd_treat(j).abund(bcd_treat(j).CXCR4vec==1);
bcd_treat(j).totabundCXCR4hi= sum(abundCXCR4hi);

indCXCR4hi = index(bcd_treat(j).CXCR4vec==1);
abundCD18hi = bcd_treat(j).abund(bcd_treat(j).CD18vec==1);
bcd_treat(j).totabundCD18hi = sum(abundCD18hi);
bcd_treat(j).totabundunk = 1-bcd_treat(j).totabundCXCR4hi-bcd_treat(j).totabundCD18hi;
indCD18hi = index(bcd_treat(j).CD18vec==1);
indneut = index(bcd_treat(j).CXCR4vec==0 & bcd_treat(j).CD18vec==0);
abundneut = bcd_treat(j).abund(indneut);

bar(indCXCR4hi, abundCXCR4hi, 'r', 'EdgeColor','none')
hold on
bar(indCD18hi, abundCD18hi, 'b', 'FaceAlpha', 0.4, 'EdgeColor','none')
bar(indneut, abundneut, 'k', 'FaceAlpha', 0.2,'EdgeColor','none')
xlabel ('barcodes')
ylabel('Abundance')
legend('CXCR4^{hi} lineages', 'CD18^{hi} lineages', 'Not known')
legend boxoff
%title('Growth rate distribution stratified by CXCR4')
set(gca,'FontSize',14,'LineWidth',1.5, 'Yscale', 'log')
title([ bcd_treat(j).rep])
end
%%
for j = 6:7
    %subplot(2,3,j)
    figure;
  index = [];
  abundCXCR4hi = [];
  indCXCR4hi = [];
  abundCD18hi = [];
  indCD18hi = [];
  indneut = [];
  abundneut = [];
 index = 1:1:length(bcd_treat(j).abund);
abundCXCR4hi = bcd_treat(j).abund(bcd_treat(j).CXCR4vec==1);
bcd_treat(j).totabundCXCR4hi= sum(abundCXCR4hi);

indCXCR4hi = index(bcd_treat(j).CXCR4vec==1);
abundCD18hi = bcd_treat(j).abund(bcd_treat(j).CD18vec==1);
bcd_treat(j).totabundCD18hi = sum(abundCD18hi);
bcd_treat(j).totabundunk = 1-bcd_treat(j).totabundCXCR4hi-bcd_treat(j).totabundCD18hi;
indCD18hi = index(bcd_treat(j).CD18vec==1);
indneut = index(bcd_treat(j).CXCR4vec==0 & bcd_treat(j).CD18vec==0);
abundneut = bcd_treat(j).abund(indneut);

bar(indCXCR4hi, (abundCXCR4hi), 'r', 'FaceAlpha', 0.7,'EdgeColor','none')
hold on
bar(indCD18hi, (abundCD18hi), 'b', 'FaceAlpha', 0.4, 'EdgeColor','none')
bar(indneut, (abundneut), 'k', 'FaceAlpha', 0.2, 'EdgeColor','none')
xlabel ('barcodes')
xlim([0 length(index)])
ylabel('Abundance')
legend('CXCR4^{hi}', 'CD18^{hi}', 'Not known')
legend boxoff
%title('Growth rate distribution stratified by CXCR4')
set(gca,'FontSize',14,'LineWidth',1.5, 'Yscale', 'log')
ylim([1e-4, 1])
title([ bcd_treat(j).rep])
end


