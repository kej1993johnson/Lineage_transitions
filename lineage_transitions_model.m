% Model of lineage resolution and phenotypic transitions in response to
% treatment
close all; clear all; clc
% This code runs simulations of the transitions between three cell types to
% reach the observed untreated proportions of 3 phenotypic clusters which we
% call here S, R, and T
%%
 tspan = [0, 528]; % time span = 2 weeks in hours
 % Initial cell numbers of S, R, and T cells where we assume 3 lineages
 % each start in these three independent cell states
 % (1) = S
 % (2) = T
 % (3) = R
 Ni(1)=1;
 Ni(2)=1;
 Ni(3)=1;
 % Initial guess at proliferation rates
 gs= 0.03;
 gt= 0.01;
 gr= 0.005;
 % Transition rates
 ksr=0.001;
 kst=0.003;
 kts=0.001;
 ktr=0.002;
 krs=0.003;
 krt= 0.002;
 params = vertcat(gs,gt,gr,ksr,kst,kts,ktr,krs,krt);
 
 save('../out/Ni.mat', 'Ni')
 save('../out/params', 'params')
 
 % Write out ODEs
 
 f = @(t, N) [ gs*N(1)+krs*N(3)+kts*N(2)-ksr*N(1)-kst*N(1); %dS/dt
            gt*N(2)+krt*N(3)+kst*N(1)-kts*N(2)-ktr*N(2); %dT/dt
            gr*N(3)+ksr*N(1)+ktr*N(2)-krs*N(3)-krt*N(3)]; %dR/dt
 
options1 = odeset('Refine',1);  
options = odeset(options1,'NonNegative',1:3);
[t,N]=ode45(f, tspan,Ni, options);

% Your output
S = N(:,1);
T = N(:,2);
R = N(:,3);
Nt = S+R+T;
fs = S./(T+R+S);
ft = T./(T+R+S);
fr = R./(T+R+S);
prop = horzcat(fs,ft,fr);
Noutgrow = N(end, :);
%% Plot figures displaying cell number over time

cells ={ 'N_{S}','N_{T}','N_{R}'};
col={ 'b','g','r'};
figure;
hold off
for i = 1:length(Ni)
plot(t,N(:,i), 'color',col{i}, 'LineWidth',2)
hold on
end
legend('N_{S}','N_{T}','N_{R}')
xlabel('time(hrs)')
ylabel('Cell number')
xlim([0, tspan(2)])
title('Number of cells in each state over time')

legend boxoff

figure;
for i = 1:length(Ni)
subplot(1,3,i)
hold off
plot(t,N(:,i), 'color', col{i}, 'LineWidth',2)
title(cells{i});
hold on
xlabel('time(hrs)')
xlim([0, tspan(2)])
ylim([0,max(N(end,:))])
ylabel('Number of cells')
end
% Plots displaying proportions of cells in time
 fracs ={ 'f_{S}','f_{T}','f_{R}'};
col={ 'b','g','r'};
figure;
hold off
for i = 1:length(Ni)
plot(t,prop(:,i), 'color',col{i}, 'LineWidth',2)
hold on
end
legend('f_{S}','f_{T}','f_{R}')
xlabel('time(hrs)')
ylabel('Proportion of cells')
xlim([0, tspan(2)])
ylim([0,1])
legend boxoff
title('Proportion of cells in each state over time: pretreatment')
 
 %% Model apply treatment as change in transition rates
  % Transition rates
 kprsr=16*ksr;
 kprst=16*kst;
 kprts=kts;
 kprtr=64*ktr;
 kprrs=krs;
 kprrt=krt;

 
 % Write out ODEs
 tspantr=[0, 12];
 f = @(t, N) [ gs*N(1)+kprrs*N(3)+kprts*N(2)-kprsr*N(1)-kprst*N(1); %dS/dt
            gt*N(2)+kprrt*N(3)+kprst*N(1)-kprts*N(2)-kprtr*N(2); %dT/dt
            gr*N(3)+kprsr*N(1)+kprtr*N(2)-kprrs*N(3)-kprrt*N(3)]; %dR/dt
 
options1 = odeset('Refine',1);  
options = odeset(options1,'NonNegative',1:3);
Nitreat= Noutgrow;
[ttr,Ntreat]=ode45(f, tspantr,Nitreat, options);
Ntot = vertcat(N,Ntreat);
ttot = vertcat(t, t(end)+ttr);


% Your output
S = Ntot(:,1);
T = Ntot(:,2);
R = Ntot(:,3);
Nt = S+R+T;
fs = S./(T+R+S);
ft = T./(T+R+S);
fr = R./(T+R+S);
prop = horzcat(fs,ft,fr);

figure;
hold off
for i = 1:length(Ni)
plot(ttot,Ntot(:,i), 'color',col{i}, 'LineWidth',2)
hold on
end
legend('N_{S}','N_{T}','N_{R}')
xlabel('time(hrs)')
ylabel('Cell number')
xlim([0, ttot(end)])
title('Number of cells in each state over time')
legend boxoff
 
figure;
for i = 1:length(Ni)
subplot(1,3,i)
hold off
plot(ttot,Ntot(:,i), 'color', col{i}, 'LineWidth',2)
title(cells{i});
hold on
xlabel('time(hrs)')
xlim([0, ttot(end)])
ylim([0,max(Ntot(end,:))])
ylabel('Number of cells')
end

figure;
hold off
for i = 1:length(Ni)
plot(ttot,prop(:,i), 'color',col{i}, 'LineWidth',2)
hold on
end
legend('f_{S}','f_{T}','f_{R}')
xlabel('time(hrs)')
ylabel('Proportion of cells')
xlim([0, ttot(end)])
ylim([0,1])
legend boxoff
title('Proportion of cells in each state over time: treat @ 408 hrs')
%%
figure;
hold off
for i = 1:length(Ni)
plot(ttot,prop(:,i), 'color',col{i}, 'LineWidth',2)
hold on
end
legend('f_{S}','f_{T}','f_{R}')
xlabel('time(hrs)')
ylabel('Proportion of cells')
xlim([528, ttot(end)])
ylim([0,1])
legend boxoff
title('Proportion of cells in each state over time: treat @ 408 hrs')
