close all; clear all; clc

[N,T] = xlsread('../data/dox_231_time.xls');
tvec = linspace (0,34,5);
colorsets = varycolor(5);
N(2:end,2:end) = N(2:end,2:end);
for i = 1:length(N)-1
    barcode(i).cellnum = N(:,1);
    barcode(i).time = N(:,i+1);
    barcode(i).well = T(i+1);
    barcode(i).color = colorsets(i, :);
    
    time = [];
    Num = [];
    time = barcode(i).time(2:end,:);
    n = length(time);
    num_params = 1;
    ind0 = find(time == 0);
    % again fit on all but t=0
    Num = barcode(i).cellnum(2:end,:); % change depend on type of processing we choose
    N0 = barcode(i).cellnum(1); % start by using what we think FACS plated
    barcode(i).Nbar = mean(barcode(i).cellnum); % average value of volumes over all measured times
    
    % First fit each curve to single exponential over entire time course
    LB = -Inf ;  % Lower Bounds
    UB = Inf; % Upper Bounds
    params0 = 0.1;% Initial Guess...
    options = optimset('TolFun',1e-12,'Tolx',1e-12,'MaxIter',1000,'Display','off','FinDiffRelStep',1e-3);
    [g, resnorm, reslsq]= lsqnonlin(@fit_singleexp, params0, LB, UB, options, Num, N0, time);
    barcode(i).g = g;
    barcode(i).residuals_g = reslsq;
    barcode(i).Rsq = 1- (sum((reslsq.^2))./(sum((barcode(i).Nbar-barcode(i).cellnum).^2)));
    barcode(i).AIC = n*log(sum(reslsq.^2)./n) +2*num_params;
    barcode(i).Nmodel = singleexpmodel(g,N0, tvec); % model fit for plotting (include t=0)
end

%%
wells={'B1','B2', 'B4', 'B5', 'B6'};
figure;
for i = 1:5
subplot(1,5,i)
plot(barcode(i).time, barcode(i).cellnum, '*', 'color', barcode(i).color, 'LineWidth',2)
hold on
plot(tvec, barcode(i).Nmodel, 'color', barcode(i).color)
xlabel ('time(days)')
ylabel('N')
title(wells(i))
xlim([0 barcode(i).time(end)])
end
% Also plot on single plot
figure;
for i = 1:5
plot(barcode(i).time, barcode(i).cellnum, '*', 'color', barcode(i).color, 'LineWidth',2)
hold on
plot(tvec, barcode(i).Nmodel, 'color', barcode(i).color)
end
xlabel ('time(days)')
ylabel('N')
title('Fitted trajectories')
ylim([ 0 2e7])

%% Make an Nmat to find the mean and variance
for i = 1:length(barcode)
    Nmat(:,i) = barcode(i).Nmodel;
end

mu_t = mean(Nmat,2);
n_2_t = mean((Nmat.^2),2);
var_t = n_2_t - ((mu_t).^2);

figure;
plot(tvec, mu_t,'b*-', 'LineWidth',2)
xlabel ('time (days)')
ylabel('mean N')
title ('mean cell number in time all wells')

figure;
plot(tvec, var_t,'r*-', 'LineWidth',2)
xlabel ('time (days)')
ylabel('variance in N')
title ('variance in cell number in time all wells')

%% Try fitting on mu_t and var_t
p = [0.7, 0.05];
% test mu and V fxns
tsamp = 0:1:tvec(end);
N0 = 180;
V0=20;
mu_n = mu_fxn(p,tsamp, N0);
V_n= V_fxn(p,tsamp,N0, V0);
figure;
plot(tsamp, mu_n, 'r')
xlabel('time (days)')
ylabel('mean N')
title('model mean')
figure;
plot(tsamp, V_n, 'b')
xlabel('time (days)')
ylabel('variance N')
title('model variance')
%%
p0 = [0.8, 0.15];
LB = [ -Inf, -Inf];
UB = [Inf, Inf];

[pfit, resnorm, res]= lsqnonlin(@fit_mu_var, p0, LB, UB, options, mu_t(2:end), var_t(2:end),tvec(2:end), N0, V0);





