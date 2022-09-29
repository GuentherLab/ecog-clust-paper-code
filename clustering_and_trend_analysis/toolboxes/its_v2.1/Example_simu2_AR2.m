%% Example of Information Dynamics for simulated bivariate system - linear estimator
clear;close all;clc;

%% Parameters
N=500; % length of simulated time series
jj=1; % index of target series
ii=2; % index of driver series
C=0.5; % Coupling parameter

%% Simulation and exact values
M=2;
p=2;
par.poles=([0.9 0.3; 0.8 0.1]); % Oscillations
par.coup=[1 2 2 C]; %in each row: "i j k c" to impose coupling from i to j at lag k with coeff c 
par.Su=[2 1]; %variance of innovation processes

%%% Theoretical VAR process
[Am,Su]=var_simulations(M,par); % parameters

[PSD,~,freq] = var_spectra(Am,Su,512,1); % spectra
PSDx=squeeze(abs(PSD(ii,ii,:)));
PSDy=squeeze(abs(PSD(jj,jj,:)));

ret = its_CElinVAR(Am,Su,jj,ii,p); % exact values information dynamics
Hy=ret.Hy; % System information
Uy_x=ret.Hy_yx; % Unexplained information
Sy=ret.Hy-ret.Hy_y; % Information Storage
Txy=ret.Hy_y-ret.Hy_yx; % Information Transfer


%% Estimation on a realization of the simulation
Un = mvnrnd(zeros(1,M),Su,N);
Yn=var_filter(Am,Un); % realization 

% embedding vector
m=p*ones(1,M);
tau=ones(1,M);
VL=its_SetLag(m,tau);

% System Information Y
eHy=its_Elin(Yn(:,jj));

% Information Storage Y
out=its_SElin(Yn,jj,VL);
eSy=out.Sy;
pval_eSy=out.p_Sy;

% Information Transfer X->Y
out=its_BTElin(Yn,ii,jj,VL);
eUy_x=out.Hy_xy;
eTxy=out.Txy;
pval_eTxy=out.p_Txy;


%% display and plot
clc;
disp('System information:');
disp(['Hy=' num2str(Hy) ' , eHy=' num2str(eHy)]);
disp(' ');
disp('Information Storage:');
disp(['Sy=' num2str(Sy) ' ; eSy=' num2str(eSy) ', p=' num2str(pval_eSy)]);
disp(' ');
disp('Information Transfer:');
disp(['Txy=' num2str(Txy) ' ; eTxy=' num2str(eTxy) ', p=' num2str(pval_eTxy)]);
disp(' ');
disp('New information:');
disp(['Ny_x=' num2str(Uy_x) ' , eNy_x=' num2str(eUy_x)]);

figure(1);
subplot(2,3,[1 2]); plot(Yn(:,ii)); zoom xon; title('Xn');
subplot(2,3,3); plot(freq,PSDx); xlim([0 0.5]); title('PSD(Xn)');
subplot(2,3,[4 5]); plot(Yn(:,jj)); zoom xon; title('Yn');
subplot(2,3,6); plot(freq,PSDy); xlim([0 0.5]); title('PSD(Yn)');


