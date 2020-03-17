% Find overlap in FACS CXCR4hi and CD18 hi barcodes with 10x high tol and
% low tol barcodes
close all; clear all; clc
%% Load in barcodes
% barcodes from FACS grouped by CXCR4 hi and low
TCXCR4hi=readtable('../out/barcodes_CXCR4_hi.csv');
TCD18hi= readtable('../out/barcodes_CD18_hi.csv');
% barcodes in 10x at TP0, with counts of high tol and low tol cells
[N,T] = xlsread('../out/hi_lo_tol.xlsx');
% make table 
lins10x = [];
for i = 1:length(N)
lins10x= vertcat(lins10x, T{i+1,1});
end
lowtolcounts = N(:,1);
hitolcounts = N(:,2);
Thilotol = table(lins10x,lowtolcounts, hitolcounts);
%% Make table
Thilotol.in10x= zeros(height(Thilotol),1);
Thilotol.inCD18hi= zeros(height(Thilotol),1);
Thilotol.inCXCR4hi= zeros(height(Thilotol),1);
Thilotol.hitol10=zeros(height(Thilotol),1);
for i = 1:length(N)
    for j = 1:height(TCD18hi)
        if char(TCD18hi.CD18hibcds(j))==Thilotol.lins10x(i,:)
            Thilotol.in10x(i)=1;
            Thilotol.inCD18hi(i) =1;
        end
    end
    for m = 1:height(TCXCR4hi)
        if char(TCXCR4hi.CXCR4hibcds(m))==Thilotol.lins10x(i,:)
            Thilotol.in10x(i)=1;
            Thilotol.inCXCR4hi(i)=1; 
        end
    end
    if Thilotol.hitolcounts(i)>Thilotol.lowtolcounts(i)
        Thilotol.hitol10x(i)=1;
    end
    %disp(i);
    
end

num_hitol10x= sum(Thilotol.hitol10x);
%%
imatch=(Thilotol.in10x==1);
Tcons = Thilotol(imatch, :);
ct_hi=0;
ct_himatch=0;
ct_lo=0;
ct_lomatch=0;
for i = 1:height(Tcons)
    if Tcons.inCXCR4hi(i)==1
        ct_hi=ct_hi+1;
        if Tcons.hitolcounts(i)>Tcons.lowtolcounts(i)
        ct_himatch=ct_himatch+1;
        end 
    end
    
    if Tcons.inCD18hi(i)==1
        ct_lo=ct_lo+1;
        if Tcons.lowtolcounts(i)>Tcons.hitolcounts(i)
        ct_lomatch=ct_lomatch+1;
        end 
    end
end    
pct_hi_match=ct_himatch/ct_hi;    
pct_lo_match=ct_lomatch/ct_lo;
            
            