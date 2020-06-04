% barcode abundace for TP1 CD18+ and TP1 CXCR4+

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
% note some barcodes are neither, these won't be categorized
end

% Make a big list of the CXCR4hi and CD18hi bcds
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
% TABLES OF CXCR4hi and CD18hi barcodes
TCXCR4= table(CXCR4hibcds);
TCD18=table(CD18hibcds);

writetable(TCXCR4, '../out/barcodes_CXCR4_hi.csv')
writetable(TCD18, '../out/barcodes_CD18_hi.csv')
%% Find which barcodes from barcode samplinf experiment are found in the 
% list of CD18 and CXCR4 positive barcodes

% This is in this order: TP0bcds, TP0 reads, rep1 bcds, rep1 reads (from
% TP50 export: rep4 bcds, rep4 reads, TP1 bcds, TP1 reads
%** Note left out the barcodes
%[N12, T12]= xlsread('../data/TP1_sorted_bcd_sampling.xlsx');
[CD18tab]= readtable('../data/CD18test.csv');
[CXCR4tab]= readtable('../data/CXCR4test.csv');


[NCD18ord, ordCD18]=sort(CD18tab.Var2, 'descend','MissingPlacement','last');
CD18bcds = CD18tab.Var1;
CD18bcds = CD18bcds(ordCD18);

[NCXCR4ord, ordCXCR4]=sort(CXCR4tab.Var2, 'descend','MissingPlacement','last');
CXCR4bcds = CXCR4tab.Var1;

%%
% TP1 not in order of abundance, so fix this

[NCD18ord, ordCD18]=sort(N12(:,1), 'descend','MissingPlacement','last');
CD18bcds = T12(2:end, 1);
T12(2:end,1) = CD18bcds(ordCD18);
N12(:,1) = NCD18ord;

[NCXCR4ord, ordCXCR4]=sort(N12(:,3), 'descend','MissingPlacement','last');
CXCR4bcds = T12(2:end, 3);
T12(2:end,3) = CDXCR4bcds(ordCXCR4);
N12(:,3) = NCXCR4ord;


%% Make a new structure for each TP/replicate
% Break this up into one structure with separate matrices
bcd_treat(1).rep='TP1 CD18';
bcd_treat(1).barcodes = CD18bcds;
bcd_treat(1).Nreads = NCD18ord;
bcd_treat(1).totreads= nansum(NCD18ord);
bcd_treat(1).abund = bcd_treat(1).Nreads./bcd_treat(1).totreads;

bcd_treat(2).rep='TP1 CXCR4';
bcd_treat(2).barcodes = CXCR4bcds;
bcd_treat(2).Nreads = NCXCR4ord;
bcd_treat(2).totreads= nansum(NCXCR4ord);
bcd_treat(2).abund = bcd_treat(2).Nreads./bcd_treat(2).totreads;



%% Test out a barcode abundance plot
figure;
bar(1:1:length(bcd_treat(1).abund), bcd_treat(1).abund)
xlabel('barcode')
ylabel('Abundance')
title('TP1 CD18+')
set(gca,'FontSize',20,'LineWidth',1.5, 'Yscale', 'Log')
%% Go through each barcode and label by CXCR4/CD18 status status
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

%% Find the median rank of CXCR4 and CD18 barcodes
for i = 1:length(bcd_treat)
    rank = [];
    CXCR4rank = [];
    CD18rank = [];
    unkrank = [];
rank = 1:1:length(bcd_treat(i).CXCR4vec);
CXCR4rank = rank(bcd_treat(i).CXCR4vec==1);
CXCR4weights = bcd_treat(i).Nreads(bcd_treat(i).CXCR4vec==1);
wCXCR4rank = (CXCR4weights).*CXCR4rank';
medCXCR4rank = median(CXCR4rank);
medwCXCR4rank = median(wCXCR4rank./bcd_treat(i).totreads);

CD18rank = rank(bcd_treat(i).CD18vec==1);
CD18weights = bcd_treat(i).Nreads(bcd_treat(i).CD18vec==1);
wCD18rank = (CD18weights).*CD18rank';
medCD18rank = median(CD18rank);
medwCD18rank = median(wCD18rank./bcd_treat(i).totreads);

unkrank = rank(bcd_treat(i).CD18vec==0 & bcd_treat(i).CXCR4vec==0);
medunkrank = median(unkrank);
unkweights = bcd_treat(i).Nreads(bcd_treat(i).CD18vec==0 & bcd_treat(i).CXCR4vec==0);
wunkrank = (unkweights).*unkrank';
medwunkrank = median(wunkrank./bcd_treat(i).totreads);

bcd_treat(i).medCXCR4rank = medCXCR4rank;
bcd_treat(i).medwCXCR4rank = medwCXCR4rank;
bcd_treat(i).medCD18rank = medCD18rank;
bcd_treat(i).medwCD18rank = medwCD18rank;
bcd_treat(i).medunkrank = medunkrank;
bcd_treat(i).medwunkrank = medwunkrank;
end
%% Try plotting barcode abundance colored by CXCR4 and CD18 status
figure;
for j = 1:length(bcd_treat)
    subplot(1,2,j)
   
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
bar(indCD18hi, abundCD18hi, 'b',  'EdgeColor','none')
bar(indneut, abundneut, 'g','EdgeColor','none')
xlabel ('barcodes')
ylabel('Abundance')
legend('TP1 CD18^{hi} lineages', 'CXCR4^{hi} lineages', 'Not known')
legend boxoff
xlim([0 100])
%title('Growth rate distribution stratified by CXCR4')
set(gca,'FontSize',14,'LineWidth',1.5)%, 'Yscale', 'log')
title([ bcd_treat(j).rep,': ', num2str(100*round(bcd_treat(j).totabundCXCR4hi, 3)),'% CXCR4^+, ', num2str(100*round(bcd_treat(j).totabundCD18hi, 3)),'% CD18^+'])
end

