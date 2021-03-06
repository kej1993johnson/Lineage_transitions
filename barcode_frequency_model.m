% Model of growth rate distributions and the effect on barcoding frequency
close all; clear all; clc
%% Start with assumption that every cell receives a unique barcode based on
N0 = 8.4e5;
q = 120e6;

[p_norepeat, Enunique, Enreps]= repeat_model(N0,q)
% Confirm this is fine....


Etot = Enunique+Enreps;
test = (Etot == N0)
prop_unique = Enunique/N0
 

%% First find apparent growth rate
% if at t=0d, sort out 5000 cells, for an initial diversity of 5000, then
% apparent growth rate from 
% N(t=24d) = N(t=0)e^(g*t)
% Nt/N0= e^(g*t)
% ln both sides
% ln(Nt)-ln(N0) = g*t
% g = (ln(Nt)-ln(N0))/(t)
N01=N0;
N11=2.4e8; % N at time of scRNA seq 
t11=24; % time of outgrowth before scRNAseq in days
g_app1 = (log(N11)-log(N01))/t11 % in units of per cell per day
tdouble1 = log(2)/g_app1
% Check what happens if we assume that this is the average growth rate:
N02=2e7;
N03 =1.5e7
N13=1.5e8;
t13=7;
g_app2 = (log(N13)-log(N03))/t13 % in units of per cell per day
tdouble2 = log(2)/g_app2
%% How does g apparent vary with sigma?
% start with your g_app1 = mu_g 
sigmavec=[];
apparent_gs = [];
sigmavec = logspace(-10, -1, 100);

for i = 1:length(sigmavec)
    sigma = sigmavec(i);
    grdistrib= [];
    grdistrib = normrnd(g_app1,sigma,N0,1);
for j = 1:length(grdistrib)
    N_g(j,i)= exp(grdistrib(j).*t11);
end

sum_Ns= sum(N_g,1);

apparent_gs(i)= (log(sum_Ns(i))-log(N0))/t11;
end

% Plot results
figure;
plot(sigmavec,apparent_gs, '*-', 'LineWidth',2)
hold on
plot(sigmavec, g_app1*ones(length(sigmavec),1), '--', 'LineWidth',2)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('standard deviation (\sigma)')
ylabel('apparent (observed) growth rate')
title('Effect of varying standard deviation on observed growth rate')
legend('apparent g', '\mu_{g}', 'Location', 'NorthWest')
legend boxoff
%% How does g apparent vary with time from barcoding to sampling?
% again start with g_app1 as your mu_g
tvecvec=[];
apparent_gs = [];
tvec = linspace(1,100, 100);
sigma = 0.05;
grdistrib = normrnd(g_app1,sigma,N0,1);
for i = 1:length(tvec)
    ti = tvec(i);
for j = 1:length(grdistrib)
    N_t(j,i)= exp(grdistrib(j).*ti);
end

sum_Ns= sum(N_t,1);

apparent_gs(i)= (log(sum_Ns(i))-log(N01))/tvec(i);
end

% Plot results
figure;
plot(tvec,apparent_gs, '*-', 'LineWidth',2)
hold on
plot(tvec, g_app1*ones(length(tvec),1), '--', 'LineWidth',2)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('Time until sampling (days)')
ylabel('apparent (observed) growth rate')
title('Effect of varying time til sampling on observed growth rate')
legend('apparent g', '\mu_{g}', 'Location', 'NorthWest')
legend boxoff

%% Create a normal distribution around this growth rate with an arbitrary
% standard deviation
stdev= 0.05;

figure;
cell_growth_rates = normrnd(g_app1,stdev,N0,1);
histogram(cell_growth_rates)
xlabel('growth rate')
xlim([0 0.6])
ylabel('frequency')
title('Normal growth rate distribution for 840,000 cells, g = 0.2525 cells/day, \sigma = 0.05, N_{final}=7.39e8')

%% Perform optimization

testparam = [ g_app1];
err_test = error_norm(testparam,stdev, N01, N11, t11)
check = sqrt(err_test)

% Fit to observed data
param0=0.22;
[g_fit1] = fminsearch(@(param)error_norm(param, stdev, N01, N11,t11),param0,optimset('TolX',1e-6, 'MaxFunEvals', 1e8));

%% Use g_fit1 to get growth rate distribution and corresponding cell number distribution

cell_growth_rates = normrnd(g_fit1,stdev,N01,1);
cell_growth_rates = sort(cell_growth_rates,'descend');

for j = 1:length(cell_growth_rates)
    N_cells(j,1)= exp(cell_growth_rates(j).*t11);
end

sum_Ncells= sum(N_cells)

Nfpop=sum_Ncells;
apparent_g = (log(Nfpop)-log(N01))/t11;
for j = 1:N01
sim_bcd_abund(j) = (N_cells(j,1)./Nfpop);
end
prob_check = sum(sim_bcd_abund)
%% Plot pdf and cell population distribution
figure;
hold off
bar(1:1:length(sim_bcd_abund), sim_bcd_abund)
%xlim([0 100])
ylabel('Relative Abundance (PDF)')
xlabel('Barcode_{i}')
%title('Simulated barcode abundance for normally distributed growth rates')
title ('Expected PDF of barcodes after t=24 days')
set(gca,'FontSize',20,'LineWidth',1.5)
N_cells_sort = sort(N_cells(:,1), 'descend');

Nchecktot = sum(N_cells_sort)

figure;
hold off
bar(1:1:length(N_cells), round(N_cells_sort,0))
set(gca, 'YScale', 'log')
%xlim([0 50])
ylabel('Number of cells in each barcode')
xlabel('Barcode_{i}')
set(gca,'FontSize',20,'LineWidth',1.5)
%title('Simulated barcode abundance for normally distributed growth rates')
title (['Total population of N_{fin} =', num2str(round(Nfpop,0)),', N_{0}=',num2str(N01),', t=', num2str(t11),' days, \mu_{g}=', num2str(round(g_fit1,2)),' cells/day']) 
%% Run estimate forward for Nsampvec at each sampling
Nsamp1 = N02;
Nvecest1 = round(Nsamp1.*sim_bcd_abund,0);
uniq_bcd_est1 = Nvecest1>=1;
N_uniq_est1 = nnz(uniq_bcd_est1)

figure;
hold off
bar(1:1:length(Nvecest1), Nvecest1)
set(gca, 'YScale', 'log')
%xlim([0 50])
ylabel('Rounded number of cells')
xlabel('Barcode_{i}')
set(gca,'FontSize',20,'LineWidth',1.5)
%title('Simulated barcode abundance for normally distributed growth rates')
title (['N_{estimated unique barcodes} =', num2str(N_uniq_est1), ', N_{samp} =', num2str(Nsamp1),', N_{0}=',num2str(N01),', N_{fin}=', num2str(N11),', t=', num2str(t11),' days, g=', num2str(round(g_fit1,2)),' cells/day']) 

% Estimate 15 million from the 20 million sample
Nsamp2 = 1.5e7;
Nvecest2 = round(((Nvecest1./Nsamp1)*Nsamp2),0);
uniq_bcd_est2 = Nvecest2>=1;
N_uniq_est2 = nnz(uniq_bcd_est2)
sumN2pop=sum(Nvecest2);

% Propagate growth forward for 7 days
t13 = 8.3;
for j = 1:length(cell_growth_rates)
    Nvecest3(j)= Nvecest2(j).*exp(cell_growth_rates(j).*t13);
end

sumN3pop= sum(Nvecest3)


apparent_g3 = (log(sumN3pop)-log(sumN2pop))/t13;
for j = 1:N01
sim_bcd_abund3(j) = (Nvecest3(j)./sumN3pop);
end
prob_check = sum(sim_bcd_abund3)

Nsamp3 = 5e6;
Nvecest3 = round(Nsamp3.*sim_bcd_abund3,0);
uniq_bcd_est3 = Nvecest3>=1;
N_uniq_est3 = nnz(uniq_bcd_est3)

figure;
hold off
bar(1:1:length(Nvecest3), Nvecest3)
set(gca, 'YScale', 'log')
%xlim([0 50])
ylabel('Rounded number of cells')
xlabel('Barcode_{i}')
set(gca,'FontSize',20,'LineWidth',1.5)
%title('Simulated barcode abundance for normally distributed growth rates')
title (['N_{estimated unique barcodes} =', num2str(N_uniq_est3), ', N_{samp} =', num2str(Nsamp3),', N_{03}=',num2str(N03),', N_{fin}=', num2str(N13),', t=', num2str(t13),' days, g=', num2str(round(g_fit1,2)),' cells/day']) 


%% First, try sampling 20 million from this...

cum_pdf = cumsum(sim_bcd_abund);
figure;
plot(1:1:length(sim_bcd_abund), cum_pdf, 'LineWidth', 2)
ylabel('CDF')
xlabel('Barcode_{i}')
set(gca,'FontSize',20,'LineWidth',1.5)
title('CDF of barcode probabilities after t=24 days')

% Sample from this distribution using the cdf
% Draw Nsamp random numbers from uniform distribution between 0 and 1

% r1 is your cell identifier
% find where this falls in the distribution
Nsampvec = zeros(N01,1);
for i = 1:Nsamp1
    rs = rand;
    ij= find(rs<cum_pdf,1, 'first');
    Nsampvec(ij) = Nsampvec(ij)+1;
end
% Check some things

Nchecksamp = sum(Nsampvec)
uniq_barcode = Nsampvec>=1;

N_uniq = nnz(uniq_barcode)

figure;
hold off
bar(1:1:length(Nsampvec), Nsampvec)
set(gca, 'YScale', 'log')
%xlim([0 50])
ylabel('Number of cells')
xlabel('Barcode_{i}')
set(gca,'FontSize',20,'LineWidth',1.5)
%title('Simulated barcode abundance for normally distributed growth rates')
title (['N_{unique barcodes} =', num2str(N_uniq), ', N_{samp} =', num2str(Nsamp1),', N_{0}=',num2str(N01),', t=', num2str(t11),' days, g=', num2str(round(g_fit1,2)),' cells/day']) 


%% Sample 1.5 million from 2 million
sim_bcd_abund2 = Nsampvec./Nsamp1;
cum_pdf2 = cumsum(sim_bcd_abund2);

Nsampvec2 = zeros(N01,1);
for i = 1:Nsamp2
    rs = rand;
    ij= find(rs<cum_pdf2,1, 'first');
    Nsampvec2(ij) = Nsampvec2(ij)+1;
end
% Check some things

Nchecksamp2 = sum(Nsampvec2)
uniq_barcode2 = Nsampvec2>=1;

N_uniq2 = nnz(uniq_barcode2)

figure;
hold off
bar(1:1:length(Nsampvec2), Nsampvec2)
set(gca, 'YScale', 'log')
%xlim([0 50])
ylabel('Number of cells')
xlabel('Barcode_{i}')
set(gca,'FontSize',20,'LineWidth',1.5)
%title('Simulated barcode abundance for normally distributed growth rates')
title (['N_{unique barcodes} =', num2str(N_uniq2),', after sampling 15 million from 20 million for viability']) 
%% Take Nsampvec and propagate forward for 8.3 days 

t13 = 8.3; 
N12 = 1.5e8% check, it should come out to near this...

for j = 1:length(cell_growth_rates)
    Nsampvec3(j)= Nsampvec2(j).*exp(cell_growth_rates(j).*t13);
end

sum_N3cells= sum(Nsampvec3)

N3fpop=sum_N3cells;
apparent_g3 = (log(N3fpop)-log(Nsamp2))/t13;
for j = 1:N01
sim_bcd_abund3(j) = (Nsampvec3(j)./N3fpop);
end
sim_bcd_abund3 = Nsampvec3./N3fpop;
cum_pdf3 = cumsum(sim_bcd_abund3);
figure;
plot(1:1:N01, cum_pdf3, 'LineWidth', 2)

figure;
hold off
bar(1:1:N01, sim_bcd_abund3)
set(gca, 'Yscale', 'log')
%%
Nsamp3 = 5e6;
Nsampvec3 = zeros(N01,1);
for i = 1:Nsamp3
    rs = rand;
    ij= find(rs<cum_pdf3,1, 'first');
    Nsampvec3(ij) = Nsampvec3(ij)+1;
end
% Check some things

Nchecksamp3 = sum(Nsampvec3)
uniq_barcode3 = Nsampvec3>=1;

N_uniq3 = nnz(uniq_barcode3)

figure;
hold off
bar(1:1:length(Nsampvec3), Nsampvec3)
set(gca, 'YScale', 'log')
%xlim([0 50])
ylabel('Number of cells')
xlabel('Barcode_{i}')
set(gca,'FontSize',20,'LineWidth',1.5)
%title('Simulated barcode abundance for normally distributed growth rates')
title (['N_{unique barcodes} =', num2str(N_uniq3),' after second outgrowth from 15 million to 150 million, sampling 5 million']) 



%% Pdf2 * Nsamp2 method of sampling
Nsamp2 = 5e6
Nvecest2 = round(Nsamp2*sim_bcd_abund2,0);
uniq_bcd_est2 = Nvecest2>=1;
N_uniq_est2 = nnz(uniq_bcd_est2)

figure;
hold off
bar(1:1:length(Nvecest2), Nvecest2)
set(gca, 'YScale', 'log')
%xlim([0 50])
ylabel('Rounded number of cells')
xlabel('Barcode_{i}')
set(gca,'FontSize',20,'LineWidth',1.5)
%title('Simulated barcode abundance for normally distributed growth rates')
title (['N_{estimated unique barcodes} =', num2str(N_uniq_est2), ', N_{samp} =', num2str(Nsamp2),', N_{0}=',num2str(N0),', t=', num2str(t11+t12),' days, g=', num2str(round(g_fit1,2)),' cells/day']) 
%% Sample from cdf of pdf2

cum_pdf2 = cumsum(sim_bcd_abund2);
figure;
plot(1:1:length(sim_bcd_abund2), cum_pdf2, 'LineWidth', 2)
ylabel('CDF')
xlabel('Barcode_{i}')
set(gca,'FontSize',20,'LineWidth',1.5)
title('CDF of barcode probabilities after t=30 days')

%% Sample from this distribution using the cdf
% Draw Nsamp random numbers from uniform distribution between 0 and 1

% r1 is your cell identifier
% find where this falls in the distribution
Nsampvec2 = zeros(N01,1);
%%
for i = 1:Nsamp2
    rs = rand;
    ij= find(rs<cum_pdf2,1, 'first');
    Nsampvec2(ij) = Nsampvec2(ij)+1;
end
%% Check some things

Nchecksamp2 = sum(Nsampvec2)
uniq_barcode2 = Nsampvec2>=1;

N_uniq2 = nnz(uniq_barcode2)

figure;
hold off
bar(1:1:length(Nsampvec2), Nsampvec2)
set(gca, 'YScale', 'log')
%xlim([0 50])
ylabel('Number of cells')
xlabel('Barcode_{i}')
set(gca,'FontSize',20,'LineWidth',1.5)
%title('Simulated barcode abundance for normally distributed growth rates')
title (['N_{unique barcodes} =', num2str(N_uniq2), ', N_{samp} =', num2str(Nsamp2),', N_{0}=',num2str(N01),', t=', num2str(t11+t12),' days, g=', num2str(round(g_fit1,2)),' cells/day']) 
%%
save('../out/Nsampvecfin.mat', 'Nsampvec3')