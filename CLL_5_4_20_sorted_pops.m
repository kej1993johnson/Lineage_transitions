% Percent Change in CXCR4 analysis (test role of CD18 cell in suppressing
% CXCR4 growth)


close all; clear all; clc;
[N, T]= xlsread('../data/4_4_20_Dilutions_for_KJ.xlsx');
names = T(2, 2:end);
nsamps = length(names);

for i = 1:nsamps
    CLLdata(i).time = N(:,1);
    CLLdata(i).sample = names(i);
    CLLdata(i).init_CXCR4 = (i-1)*10;
    CLLdata(i).init_CD18 = 100-CLLdata(i).init_CXCR4;
    CLLdata(i).CXCR4pct = N(:, i+1);
    CLLdata(i).pct_change = CLLdata(i).CXCR4pct(end)-CLLdata(i).CXCR4pct(1);
end
pct_change_loCXCR4 = [];
pct_change_hiCXCR4 = [];
figure;
for i = 2:nsamps-1
    plot(CLLdata(i).init_CXCR4, CLLdata(i).pct_change, '*')
    hold on
    if CLLdata(i).init_CXCR4<=20
        pct_change_loCXCR4 = vertcat(pct_change_loCXCR4,CLLdata(i).pct_change)
    end
    if CLLdata(i).init_CXCR4>=30
        pct_change_hiCXCR4 = vertcat(pct_change_hiCXCR4, CLLdata(i).pct_change)
    end
end
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('% CXCR4 initially')
ylabel('% increase in CXCR4')

% Are the two groups statistically significantly different?
[h,p1] = ttest2(pct_change_loCXCR4, pct_change_hiCXCR4)
