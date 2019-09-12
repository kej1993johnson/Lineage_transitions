% Barcode sampling-- code to determine how accurate of a snapshot we get of
% a barcode distribution as a function of sampling percent/size. Want to
% minimize the amount we sample while maximizing the relevant information
% obtained from the sample.
close all; clear all; clc
%% Load in OG and TPO barcode frequency distribution
s = tdfread('../data/FMV1_OG_TP0distrib.tsv', ' ', ' ');
labels = s.FM0x2D1; % this denotes your sample (i.e FM-1 to FM-7, FM-V-1-7 & OG & TP0)
abundance = s.x14;
barcodes = s.AAAAACAGAAAACGAAGACT;
%% Make a barcode structure that holds the barcode frequency distribution for each sample.
bcd = struct( 'sample', strings(1), 'n_cells', zeros(1,1),'n_unique_barcodes', zeros(1,1));
indOG = [];
indTP0 = [];
% Find the indices that correspond to your OG and TP0 samples
for i = 1:length(labels)
    if ismember(labels(i,:), 'OP     ')
        indOG = vertcat(indOG,i);
    end
    if ismember(labels(i,:), 'TP0    ')
        indTP0 = vertcat(indTP0,i);
    end
end

%% Make a vector that holds the barcode frequency distribution
distribOG = abundance(indOG);
distribTP0 = abundance(indTP0);
[distribOGord,iOG] = sort(distribOG,'descend');
[distribTP0ord,iTP0] = sort(distribTP0, 'descend');

figure;
bar(1:1:length(indOG), distribOGord)
xlabel('barcode')
ylabel('abundance')
title('OG barcode distribution')

figure;
bar(1:1:length(indTP0), distribTP0ord)
xlabel('barcode')
ylabel('abundance')
title('TP0 barcode distribution')

bcdsOG = barcodes(indOG, :);
bcdsTP0 = barcodes(indTP0, :);

% order these
bcdsOGord=bcdsOG(iOG);
bcdsTP0ord = bcdsTP0(iTP0);
% Put into your structure
bcd(1).sample = 'OG';
bcd(2).sample = 'TP0';
bcd(1).n_unique_barcodes = length(indOG);
bcd(2).n_unique_barcodes = length(indTP0);
bcd(1).abund = distribOGord;
bcd(2).abund = distribTP0ord;


for i = 1:2
    bcd(i).n_cells = sum(bcd(i).abund);
end
%% Now that you have these, simulate sampling from them
for i = 1:2
    bcd(i).pdf = bcd(i).abund./bcd(i).n_cells;
    bcd(i).cdf = cumsum(bcd(i).pdf);
end
 
%% Simulate 1 rouns of sampling 1% of the population
for i = 1:2
    bcd(i).Nsamp = round(0.01*bcd(i).n_cells,0);
    Nsampvec = zeros(bcd(i).n_unique_barcodes,1); % this is your vector to keep track of number of cells in each barcode
    for j = 1:bcd(i).Nsamp % for each sample you draw
        rs = rand;
        ij = find(rs< bcd(i).cdf, 1, 'first');
        Nsampvec(ij) = Nsampvec(ij)+1; % keep track of cells with each barcode
    end
    bcd(i).sampabund = Nsampvec;
    bcd(i).samppdf = Nsampvec./bcd(i).Nsamp; % this is your sample abundance
    bcd(i).n_unique_samp = nnz(Nsampvec>=1);
    bcd(i).CCCsamp = f_CCC([bcd(i).pdf, bcd(i).samppdf]);
    bcd(i).prop_uniq = bcd(i).n_unique_samp/bcd(i).n_unique_barcodes;
    
end
%% Plot the results
figure;
bar(1:1:length(indOG), bcd(1).pdf)
hold on
bar(1:1:length(indOG),bcd(1).samppdf)
legend('true pdf', 'sampled pdf')
xlabel('barcode')
ylabel('abundance')
title(['OG barcode distribution comparison, CCC=', num2str(bcd(1).CCCsamp)])

figure;
bar(1:1:length(indTP0), bcd(2).pdf)
hold on
bar(1:1:length(indTP0),bcd(2).samppdf)
legend('true pdf', 'sampled pdf')
xlabel('barcode')
ylabel('abundance')
title(['TP0 barcode distribution comparison, CCC=', num2str(bcd(2).CCCsamp)])
%% Simulate sampling many times at a certain percent of the population
% Calculate average percent error at each unique barcode index.

nruns = 100;

for i = 1:2
    Nsamp= round(0.01*bcd(i).n_cells,0);
     % this is your vector to keep track of number of cells in each barcode
    samppdf = [];
    pct_err = [];
    for s = 1:nruns
       Nsampvec = zeros(bcd(i).n_unique_barcodes,1);
    for j = 1:bcd(i).Nsamp % for each sample you draw
        rs = rand;
        ij = find(rs< bcd(i).cdf, 1, 'first');
        Nsampvec(ij) = Nsampvec(ij)+1; % keep track of cells with each barcode
    end
    samppdf(:,s) = Nsampvec./bcd(i).Nsamp; % this is your sample abundance  
    end
bcd(i).sampmat = samppdf; % this saves the proportion of cells in that lineage for each sample

% Find percent error at each barcode
    for s = 1:nruns
    pct_err(:,s)=abs(bcd(i).pdf-samppdf(:,s))./bcd(i).pdf;
    end
    
bcd(i).mean_pct_err = mean(pct_err,2);
end

%%
figure;
plot(1:1:length(indOG),100*bcd(1).mean_pct_err, '.')
xlabel('Barcode')
ylabel('Average error')
ylim([ 0 100])
title('Average error at each barcode OG for 1% sampling')
set(gca,'FontSize',20,'LineWidth',1.5)

figure;
plot(1:1:length(indTP0),100*bcd(2).mean_pct_err, '.')
xlabel('Barcode')
ylabel('Average % error')
ylim([ 0 100])
title('Average error at each barcode TP0 for 1% sampling')
set(gca,'FontSize',20,'LineWidth',1.5)

figure;
plot(1:1:length(indOG),100*bcd(1).mean_pct_err, '.')
hold on
plot(1:1:length(indTP0),100*bcd(2).mean_pct_err, '.')
xlabel('Barcode')
ylabel('Average % error')
legend ('OG', 'TP0')
legend boxoff
%ylim([ 0 100])
title('Average error in proportion at each barcode for 1% sampling')
set(gca,'FontSize',20,'LineWidth',1.5)

%%
figure;
plot(bcd(1).pdf,100*bcd(1).mean_pct_err)
hold on
plot(bcd(2).pdf, 100*bcd(2).mean_pct_err)
xlabel('Proportion of cells with that barcode in true population')
ylabel('Average  % error')
legend('OG', 'TP0')
ylim([ 0 100])
title('Error versus barcode abundance')
set(gca,'FontSize',20,'LineWidth',1.5)
%% Simulate sampling many times at varying percent of the population
% Calculate average percent error at each unique barcode index.
sampvec = [ .01, .02, .03, .04, .05,.06,.07, .08, .09, .1];
nruns = 10;
for i = 1:2
    mean_err_mat = [];
for m = 1:length(sampvec)

    Nsamp= round(sampvec(m)*bcd(i).n_cells,0);
    bcd(i).Nsamps(m) = Nsamp;
     % this is your vector to keep track of number of cells in each barcode
    samppdf = [];
    pct_err = [];
    
    for s = 1:nruns
       Nsampvec = zeros(bcd(i).n_unique_barcodes,1);
    for j = 1:bcd(i).Nsamps(m) % for each sample you draw
        rs = rand;
        ij = find(rs< bcd(i).cdf, 1, 'first');
        Nsampvec(ij) = Nsampvec(ij)+1; % keep track of cells with each barcode
    end

    samppdf(:,s) = Nsampvec./bcd(i).Nsamps(m); % this is your sample abundance  
    end


    bcd(i).sampmat = samppdf; % this saves the proportion of cells in that lineage for each sample

% Find percent error at each barcode
    for s = 1:nruns
    pct_err(:,s)=abs(bcd(i).pdf-samppdf(:,s))./bcd(i).pdf;
    end
    
    mean_pct_err = mean(pct_err,2);
    mean_err_mat(:,m) = mean_pct_err;
    
end
bcd(i).mean_err_vs_pct = mean_err_mat;

end
%% For each sampling percentage, find the index at which the average error goes above 50%
for i = 1:2
    mean_err_vs_pct = bcd(i).mean_err_vs_pct;
    vec = [];
    num_id_barcodes = [];
for m = 1:length(sampvec)
    vec = mean_err_vs_pct(:,m);
    num_id_barcodes(m,1) = find(vec>0.5, 1, 'first');
end
bcd(i).num_id_barcodes = num_id_barcodes;
end

%% Find the smoothed mean err vs pct
for i = 1:2
    mean_err_vs_pct = bcd(i).mean_err_vs_pct;
    vec = [];
    smoothedmeanerr = [];
    for m = 1:length(sampvec)
        vec = mean_err_vs_pct(:,m);
        smoothedmeanerr(:,m) = smoothdata(vec, 'gaussian', 100);
    end
    bcd(i).smooth_mean_err = smoothedmeanerr;
end

figure;
subplot(1,2,1)
for m = 1:length(sampvec)
plot(1:1:length(indOG),100*bcd(1).smooth_mean_err(:,m), '.')
hold on
end
legend ( '0.01%','1%','3%',  '5%',  '8%', 'Location', 'Northwest')
legend boxoff
xlabel('Barcode')
ylabel('Smoothed average % error')
ylim([ 0 500])
% xlim([ 0 2200])
title('OG')
set(gca,'FontSize',20,'LineWidth',1.5)
subplot(1,2,2)
for m = 1:length(sampvec)
plot(1:1:length(indTP0),100*bcd(2).smooth_mean_err(:,m), '.')
hold on
end
legend ( '0.01%','1%','3%',  '5%', '8%', 'Location', 'Northwest')
legend boxoff
ylim([ 0 500])
% xlim([ 0 2200])
xlabel('Barcode')
ylabel('Smoothed average % error')
title('TP0')
set(gca,'FontSize',20,'LineWidth',1.5)


%% Plot results
figure;
subplot(1,2,1)
for m = 1:length(sampvec)
plot(1:1:length(indOG),100*bcd(1).mean_err_vs_pct(:,m), '.')
hold on
end
legend ( '0.01%','1%','3%',  '5%', '8%', 'Location', 'Northwest')
legend boxoff
xlabel('Barcode')
ylabel('Average % error')
ylim([ 0 100])
xlim([ 0 2200])
title('OG')
set(gca,'FontSize',20,'LineWidth',1.5)
subplot(1,2,2)
for m = 1:length(sampvec)
plot(1:1:length(indTP0),100*bcd(2).mean_err_vs_pct(:,m), '.')
hold on
end
legend ( '0.01%','1%','3%',  '5%', '8%', 'Location', 'Northwest')
legend boxoff
% ylim([ 0 100])
% xlim([ 0 2200])
xlabel('Barcode')
ylabel('Average % error')
title('TP0')
set(gca,'FontSize',20,'LineWidth',1.5)

figure;
for i = 1:2
plot(100*sampvec, bcd(i). num_id_barcodes, 'LineWidth', 4)
hold on
end
legend('OG', 'TP0')
legend boxoff
xlabel('Percent of population sampled')
ylabel('Number of barcodes whose abundance error < 50%')
title('Number of reliable barcodes verus sampling percent')
set(gca,'FontSize',20,'LineWidth',1.5)

%% Simulate sampling 1% 3 times vs sampling 3%
% Calculate average percent error at each unique barcode index.
sampvec = [.01, .01, .01, .03];
nruns = 1;
for i = 1:2
    mean_err_mat = [];
for m = 1:length(sampvec)

    Nsampm= round(sampvec(m)*bcd(i).n_cells,0);
   
     % this is your vector to keep track of number of cells in each barcode
    samppdf = [];
    pct_err = [];
    
    for s = 1:nruns
       Nsampvec = zeros(bcd(i).n_unique_barcodes,1);
    for j = 1:Nsampm % for each sample you draw
        rs = rand;
        ij = find(rs< bcd(i).cdf, 1, 'first');
        Nsampvec(ij) = Nsampvec(ij)+1; % keep track of cells with each barcode
    end

    samppdf(:,s) = Nsampvec./Nsampm; % this is your sample abundance  
    end
    mean_pdf = mean(samppdf,2);

    bcd(i).pdfs13(:,m) = mean_pdf; % this saves the proportion of cells in that lineage for each sample

% Find percent error at each barcode
    
    pct_err=abs(bcd(i).pdf-mean_pdf)./bcd(i).pdf;
 
    bcd(i).err13(:,m) = pct_err;
    
end


end
%% Plot the results of the 1% sampling
figure;
    subplot(1,2,1)
    semilogy(1:1:length(bcd(1).pdf), bcd(1).pdf, 'k.', 'LineWidth',3)
    hold on
    for m = 1:length(sampvec)-1
    semilogy(1:1:length(bcd(1).pdf), bcd(1).pdfs13(:,m), 'g.')
    end
    for m = length(sampvec)
    semilogy(1:1:length(bcd(1).pdf), bcd(1).pdfs13(:,m), 'b.')
    end
    semilogy(1:1:length(bcd(1).pdf), bcd(1).pdf, 'k.', 'LineWidth',3)
xlabel('Barcode')
ylabel('Proportion of cells in each barcode')
legend('True pdf', '1%', '1%', '1%', '3%')
legend boxoff
title('OG pdfs for different sample sizes')
set(gca,'FontSize',20,'LineWidth',1.5)

    subplot(1,2,2)
    semilogy(1:1:length(bcd(2).pdf), bcd(2).pdf, 'k.', 'LineWidth', 3)
    hold on
    for m = 1:length(sampvec)-1
    semilogy(1:1:length(bcd(2).pdf), bcd(2).pdfs13(:,m), 'g.')
    end
    for m = length(sampvec)
    semilogy(1:1:length(bcd(2).pdf), bcd(2).pdfs13(:,m), 'b.')
    end
    semilogy(1:1:length(bcd(2).pdf), bcd(2).pdf, 'k.', 'LineWidth', 3)
xlabel('Barcode')
ylabel('Proportion of cells in each barcode')
legend('True pdf', '1%', '1%', '1%', '3%')
legend boxoff
title('TP0 pdfs for different sample sizes')
set(gca,'FontSize',20,'LineWidth',1.5)
%% Find the standard deviation and std dev of log of the 3 1% samples
for i = 1:2
    bcd(i).pdfstdev = std(bcd(i).pdfs13(:,1:3),0,2);
    bcd(i).pdflogstdev= std(log(bcd(i).pdfs13(:,1:3)),0,2);
    bcd(i).pdfcov = (std(bcd(i).pdfs13(:,1:3),0,2))./mean(bcd(i).pdfs13(:,1:3),2);
    bcd(i).pdfmean1s = mean(bcd(i).pdfs13(:,1:3),2);
end

figure;
plot (1:1:length(bcd(1).pdf), bcd(1).pdfcov, '.')
hold on
xlabel('barcode')
ylabel('COV')
set(gca,'FontSize',20,'LineWidth',1.5)
plot (1:1:length(bcd(2).pdf), bcd(2).pdfcov, '.')
legend('OG', 'TP0')
legend boxoff
title('COV in three 1% samples at each barcode')
set(gca,'FontSize',20,'LineWidth',1.5)

figure;
plot (1:1:length(bcd(1).pdf), bcd(1).pdfstdev, '.')
hold on
xlabel('barcode')
ylabel('Standard Deviation')
set(gca,'FontSize',20,'LineWidth',1.5)
plot (1:1:length(bcd(2).pdf), bcd(2).pdfstdev, '.')
legend('OG', 'TP0')
legend boxoff
title('Standard Deviation in 3 1% samples at each barcode')
set(gca,'FontSize',20,'LineWidth',1.5)
figure;
plot (1:1:length(bcd(1).pdf), bcd(1).pdflogstdev, '.')
hold on
xlabel('barcode')
ylabel('Standard deviation of log of abundance')
set(gca,'FontSize',20,'LineWidth',1.5)
plot (1:1:length(bcd(2).pdf), bcd(2).pdflogstdev, '.')
legend('OG', 'TP0')
legend boxoff
title('Standard Deviation of log of abundance in 3 1% samples at each barcode')
set(gca,'FontSize',20,'LineWidth',1.5)

figure;
%subplot(1,2,1)
plot3(bcd(1).pdflogstdev, log(mean(bcd(1).pdfs13(:,1:3),2)),100*bcd(1).err13(:,4), '.')
hold on
plot3(bcd(2).pdflogstdev, log(mean(bcd(2).pdfs13(:,1:3),2)),100*bcd(2).err13(:,4), '.')
zlabel('Percent error')
ylabel('log(mean abundance)')
xlabel('Stdev(log(abundance))')
set(gca,'FontSize',20,'LineWidth',1.5)
zlim([ 0 100])
title('Expected error as a function of abundance and spread')
legend('OG', 'TP0')
legend boxoff
%%
figure;
%subplot(1,2,1)
plot3(bcd(1).pdfcov, log(bcd(1).pdfmean1s),100*bcd(1).err13(:,4), '.')
hold on
plot3(bcd(2).pdfcov, log(bcd(2).pdfmean1s),100*bcd(2).err13(:,4), '.')
zlabel('Percent error')
ylabel('Log(mean abundance)')
xlabel('COV')
set(gca,'FontSize',20,'LineWidth',1.5)
zlim([ 0 100])
title('Expected error as a function of abundance and spread')
legend('OG', 'TP0')
legend boxoff

figure;
plot(log(bcd(1).pdfmean1s),bcd(1).pdfcov, '.')
hold on
plot(log(bcd(2).pdfmean1s),bcd(2).pdfcov, '.')
xlabel('log(mean abundance)')
ylabel('COV')
set(gca,'FontSize',20,'LineWidth',1.5)
title('COV versus abundance')

figure;
%subplot(1,2,1)
plot(bcd(1).pdfcov, 100*bcd(1).err13(:,4), '.')
hold on
plot(bcd(2).pdfcov, 100*bcd(2).err13(:,4), '.')
ylabel('Percent error')
xlabel('COV')
set(gca,'FontSize',20,'LineWidth',1.5)
zlim([ 0 100])
title('Expected error as a function of COV')
legend('OG', 'TP0')
legend boxoff
%%
figure;
subplot(1,2,1)
for m = 1:length(sampvec)-1
plot(1:1:length(bcd(1).pdf),100*bcd(1).err13(:,m), 'g.')
hold on
end
for m = length(sampvec)
plot(1:1:length(bcd(1).pdf),100*bcd(1).err13(:,m), 'b.')
hold on
end
legend ( '1%', '1%', '1%', '3%', 'Location', 'Northwest')

legend boxoff
xlabel('Barcode')
ylabel('Average % error')
title('OG average error for different sample sizes')
ylim([ 0 100])
% xlim([ 0 2200])
set(gca,'FontSize',20,'LineWidth',1.5)

subplot(1,2,2)
for m = 1:length(sampvec)-1
plot(1:1:length(bcd(2).pdf),100*bcd(2).err13(:,m), 'g.')
hold on
end
for m = length(sampvec)
plot(1:1:length(bcd(2).pdf),100*bcd(2).err13(:,m), 'b.')
hold on
end
legend ( '1%', '1%', '1%', '3%', 'Location', 'Northwest')
legend boxoff
xlabel('Barcode')
ylabel('Average % error')
title('TP0 average error for different sample sizes')
ylim([ 0 100])
% xlim([ 0 2200])
set(gca,'FontSize',20,'LineWidth',1.5)
%%
figure;
    subplot(1,2,1)
    semilogy(1:1:length(bcd(1).pdf), bcd(1).pdf, 'k.', 'LineWidth',3)
    hold on

    semilogy(1:1:length(bcd(1).pdf), mean(bcd(1).pdfs13(:,1:3),2), 'g.')

    semilogy(1:1:length(bcd(1).pdf), bcd(1).pdfs13(:,m), 'b.')
 
    semilogy(1:1:length(bcd(1).pdf), bcd(1).pdf, 'k.', 'LineWidth',3)
xlabel('Barcode')
ylabel('Proportion of cells in each barcode')
legend('True pdf', 'avg 3 1%', '3%')
legend boxoff
title('OG pdfs for different sample sizes')
set(gca,'FontSize',20,'LineWidth',1.5)

    subplot(1,2,2)
    semilogy(1:1:length(bcd(2).pdf), bcd(2).pdf, 'k.', 'LineWidth', 3)
    hold on

    semilogy(1:1:length(bcd(2).pdf), mean(bcd(2).pdfs13(:,1:3),2), 'g.')

    semilogy(1:1:length(bcd(2).pdf), bcd(2).pdfs13(:,m), 'b.')

    semilogy(1:1:length(bcd(2).pdf), bcd(2).pdf, 'k.', 'LineWidth', 3)
xlabel('Barcode')
ylabel('Proportion of cells in each barcode')
legend('True pdf', 'avg 3 1%', '3%')
legend boxoff
title('TP0 pdfs for different sample sizes')
set(gca,'FontSize',20,'LineWidth',1.5)