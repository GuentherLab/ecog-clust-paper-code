%% Example of Information Dynamics for simulated monovariate system - linear estimator
clear;close all;clc;

%% Parameters
N=500; % length of simulated time series
r=0.8; % strength of stochastic oscillation
f=0.3; % frequency of stochastic oscillation

%% Simulation and exact values
M=1; %n. of time series
p=2; % maximum lag
par.poles=([r f]); % Oscillation
par.coup=[]; %in each row: "i j k c" to impose coupling from i to j at lag k with coeff c 
par.Su=[1]; %variance of innovation processes

%%% Theoretical VAR process
[Am,Su]=var_simulations(M,par); % parameters

[PSD,~,freq] = var_spectra(Am,Su,512,1); % spectra
PSDy=squeeze(real(PSD));

ret = its_CElinVAR1(Am,Su,p); % exact values information dynamics
Hy=ret.Hy; % System information
Ny=ret.Hy_y; % New information
Sy=Hy-Ny; % Information Storage


%% Estimation on a realization of the simulation
Un = mvnrnd(zeros(1,M),Su,N);
Yn=var_filter(Am,Un); % realization 

% embedding vector
m=p*ones(1,M);
tau=ones(1,M);
VL=its_SetLag(m,tau);

% System Information
eHy=its_Elin(Yn);

% Information Storage
out=its_SElin(Yn,1,VL);
eNy=out.Hy_y;
eSy=out.Sy;
pval_eSy=out.p_Sy;

%% display and plot
clc;
disp('System information:');
disp(['Hy=' num2str(Hy) ' ; eHy=' num2str(eHy)]);
disp(' ');
disp('New information:');
disp(['Ny=' num2str(Ny) ' ; eNy=' num2str(eNy)]);
disp(' ');
disp('Information Storage:');
disp(['Sy=' num2str(Sy) ' ; eSy=' num2str(eSy) ', p=' num2str(pval_eSy)]);

figure(1);
subplot(1,3,[1 2]); plot(Yn); zoom xon; title('Yn');
subplot(1,3,3); plot(freq,PSDy); axis([0 0.5 0 1.1*max(PSDy)]); title('PSD(Yn)');


