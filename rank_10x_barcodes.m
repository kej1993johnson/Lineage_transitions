%Find median rank of the 10x barcodes separated by CXCR4, CD18 status from
%each of the three timepoitns in the 10x data

close all; clear all; clc;

%% We will use the old barcode sampling data to identify which lineages are
% CXCR4+ and which are CD8+

[N, T]= xlsread('../data/Barcode_sampling_pre_treat.xlsx');
N(isnan(N))=0;
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
%% Find which barcodes from 10x lineages (TP0, TP1, TP2) 
% list of CD18 and CXCR4 positive barcodes
[N1, T1] = xlsread('../data/lineage_timepoint_count_table.xlsx');
N1(isnan(N1))=0;
%%
[N0ord, ordN0]=sort(N1(:,1), 'descend');
bcd10x = T1(2:end,1);
bcdTP0 = bcd10x(ordN0);
bcdTP0nnz = bcdTP0(N0ord~=0);
N0ordnnz = N0ord(N0ord~=0);

[N50ord, ordN50]=sort(N1(:,2), 'descend');
bcdTP50 = bcd10x(ordN50);
bcdTP50nnz = bcdTP50(N50ord~=0);
N50ordnnz = N50ord(N50ord~=0);

[N1ord, ordN1]=sort(N1(:,3), 'descend');
bcdTP1 = bcd10x(ordN1);
bcdTP1nnz = bcdTP1(N1ord~=0);
N1ordnnz = N1ord(N1ord~=0);



%% Make a new structure for each TP/replicate
% Break this up into one structure with separate matrices
bcd_treat(1).rep='TP0 10x';
bcd_treat(1).barcodes = bcdTP0;
bcd_treat(1).bcdsnnz = bcdTP0nnz;
bcd_treat(1).Nreads = N0ord(:,1);
bcd_treat(1).Nreadsnnz = N0ordnnz;
bcd_treat(1).totreads= sum(bcd_treat(1).Nreads);
bcd_treat(1).abund = bcd_treat(1).Nreads./bcd_treat(1).totreads;
bcd_treat(2).rep = 'TP50 10x';
bcd_treat(2).barcodes = bcdTP50;
bcd_treat(2).bcdsnnz = bcdTP50nnz;
bcd_treat(2).Nreads = N50ord(:,1);
bcd_treat(2).Nreadsnnz = N50ordnnz;
bcd_treat(2).totreads= sum(bcd_treat(2).Nreads);
bcd_treat(2).abund = bcd_treat(2).Nreads./bcd_treat(2).totreads;
bcd_treat(3).rep = 'TP1 10x';
bcd_treat(3).barcodes = bcdTP1;
bcd_treat(3).bcdsnnz = bcdTP1nnz;
bcd_treat(3).Nreads = N1ord(:,1);
bcd_treat(3).Nreadsnnz = N1ordnnz;
bcd_treat(3).totreads= sum(bcd_treat(3).Nreads);
bcd_treat(3).abund = bcd_treat(3).Nreads./bcd_treat(3).totreads;


%% Go through each barcode and label by CXCR4/CD18 status status
for j = 1:length(bcd_treat)
    % default to 0 in status
    CXCR4vec = [];
    CD18vec = [];
    CXCR4vec = zeros(length(bcd_treat(j).bcdsnnz),1);
    CD18vec = zeros(length(bcd_treat(j).bcdsnnz),1);
    % if barcode matches one of the list of CXCR4/CD18 barcodes, label
    for i = 1:length(bcd_treat(j).bcdsnnz)
        if any(strcmp(TCXCR4.CXCR4hibcds,bcd_treat(j).bcdsnnz(i)))
            CXCR4vec(i)=1;
        end
        if any(strcmp(TCD18.CD18hibcds,bcd_treat(j).bcdsnnz(i)))
            CD18vec(i)=1;
        end
    end
    bcd_treat(j).CXCR4vec = CXCR4vec;
    bcd_treat(j).CD18vec = CD18vec;
end

%% Try plotting barcode abundance colored by CXCR4 and CD18 status
figure;
for j = 1:length(bcd_treat)
    subplot(2,3,j)
   
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
legend('CXCR4^{hi} lineages', 'CD18^{hi} lineages', 'Not known')
legend boxoff
xlim([0 100])
%title('Growth rate distribution stratified by CXCR4')
set(gca,'FontSize',14,'LineWidth',1.5)%, 'Yscale', 'log')
title([ bcd_treat(j).rep,' ', num2str(100*round(bcd_treat(j).totabundCXCR4hi, 3)),'% CXCR4^+'])
end



%% Using weighted rank
for i = 1:length(bcd_treat)
    rank = [];
    CXCR4rank = [];
    CD18rank = [];
    unkrank = [];
rank = 1:1:length(bcd_treat(i).CXCR4vec);
CXCR4rank = rank(bcd_treat(i).CXCR4vec==1);
CXCR4weights = bcd_treat(i).Nreadsnnz(bcd_treat(i).CXCR4vec==1);
wCXCR4rank = (CXCR4weights).*CXCR4rank';
medCXCR4rank = median(CXCR4rank);
medwCXCR4rank = median(wCXCR4rank./bcd_treat(i).totreads);

CD18rank = rank(bcd_treat(i).CD18vec==1);
CD18weights = bcd_treat(i).Nreadsnnz(bcd_treat(i).CD18vec==1);
wCD18rank = (CD18weights).*CD18rank';
medCD18rank = median(CD18rank);
medwCD18rank = median(wCD18rank./bcd_treat(i).totreads);

unkrank = rank(bcd_treat(i).CD18vec==0 & bcd_treat(i).CXCR4vec==0);
medunkrank = median(unkrank);
unkweights = bcd_treat(i).Nreadsnnz(bcd_treat(i).CD18vec==0 & bcd_treat(i).CXCR4vec==0);
wunkrank = (unkweights).*unkrank';
medwunkrank = median(wunkrank./bcd_treat(i).totreads);

bcd_treat(i).medCXCR4rank = medCXCR4rank;
bcd_treat(i).medwCXCR4rank = medwCXCR4rank;
bcd_treat(i).medCD18rank = medCD18rank;
bcd_treat(i).medwCD18rank = medwCD18rank;
bcd_treat(i).medunkrank = medunkrank;
bcd_treat(i).medwunkrank = medwunkrank;
end




%% Find the median rank of CXCR4 and CD18 cells, ordered by barcodes

for i = 1%:length(bcd_treat)
    rank = [];
    CXCR4rank = [];
    CD18rank = [];
    unkrank = [];
rank = 1:1:(bcd_treat(i).totreads);
ranklist = 0;
ranklistCXCR4 = 0;
ranklistCD18 = 0;
ranklistunk = 0;
    for j = 1:length(bcd_treat(i).bcdsnnz)
        newvals = [];
        newvals = 1:1:bcd_treat(i).Nreadsnnz(j);
        ranklist = horzcat(ranklist, ranklist(end)+newvals);
        if bcd_treat(i).CXCR4vec(j)==1
            ranklistCXCR4 = horzcat(ranklistCXCR4, ranklist(end)+newvals);
        end
        if bcd_treat(i).CD18vec(j)==1
            ranklistCD18= horzcat(ranklistCD18, ranklist(end)+newvals);
        end
        if bcd_treat(i).CXCR4vec(j)==0 && bcd_treat(i).CD18vec(j)==0
            ranklistunk = horzcat(ranklistunk, ranklist(end)+newvals);
        end
    end
        
 medCXCR4rank = median(ranklistCXCR4(2:end));
 
 medCD18rank = median(ranklistCD18(2:end));

 medunkrank = median(ranklistunk(2:end));
bcd_treat(i).medCXCR4rank = medCXCR4rank;
bcd_treat(i).medCD18rank = medCD18rank;
bcd_treat(i).medunkrank = medunkrank;
end
%%
figure;

plot(ranklistCD18(2:end), ones(length(ranklistCD18(2:end)),1), '.','color', 'b')
hold on
plot(ranklistCXCR4(2:end), ones(length(ranklistCXCR4(2:end)),1), '.','color', 'r')
plot(ranklistunk(2:end), ones(length(ranklistunk(2:end)),1), '.','color', 'g')
plot([bcd_treat(1).medCXCR4rank bcd_treat(1).medCXCR4rank], [1 2], '-','color', 'r')
plot([bcd_treat(1).medCD18rank bcd_treat(1).medCD18rank], [1 2], '-','color', 'b')
plot([bcd_treat(1).medunkrank bcd_treat(1).medunkrank], [1 2], '-','color', 'g')
set(gca,'FontSize',14,'LineWidth',1.5)
xlabel('cells')
ylim([ 0.5 2.5])
legend('CD18^{hi}', 'CXCR4^{hi}', 'unkown')
title('Cells ordered by barcode abundance')

%% figure;
figure;
for j = 1%:length(bcd_treat)
 
   
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
legend('CXCR4^{hi} lineages', 'CD18^{hi} lineages', 'Not known')
legend boxoff
%xlim([0 100])
%title('Growth rate distribution stratified by CXCR4')
set(gca,'FontSize',14,'LineWidth',1.5)%, 'Yscale', 'log')
title([ bcd_treat(j).rep,' ', num2str(100*round(bcd_treat(j).totabundCXCR4hi, 3)),'% CXCR4^+'])
end
