%% Example of Information Dynamics for simulated 3-variate system - linear estimator
clear;close all;clc;

%% Parameters
N=500; % length of simulated time series
jj=2; % index of target series Y (RR)
ii=1; % index of driver series X (resp)
C=0.5; % Coupling parameter

%% Simulation and exact values
M=3;
p=2;
kk=setdiff(1:M,[ii jj]); % index of driver series Z (sap)
par.poles=([0.9 0.3; 0.8 0.1; 0.8 0.1]); % Oscillations
% par.poles=([0.9 0.3; 0.8 0.1; 0.8 0.1]); % Oscillations
par.coup=[1 2 1 C; 1 3 2 1-C; 3 2 1 0.5; 2 3 2 0.1]; %in each row: "i j k c" to impose coupling from i to j at lag k with coeff c 
par.Su=[5 1 1]; %variance of innovation processes

%%% Theoretical VAR process
[Am,Su]=var_simulations(M,par); % parameters

[PSD,~,freq] = var_spectra(Am,Su,512,1); % spectra
PSDx=squeeze(abs(PSD(ii,ii,:)));
PSDy=squeeze(abs(PSD(jj,jj,:)));
PSDz=squeeze(abs(PSD(kk,kk,:)));

ret = its_CElinVAR(Am,Su,jj,ii,p); % exact values information dynamics
Hy=ret.Hy; % System information Y
Uy_xz=ret.Hy_yzx; % Unexplained information Y|X,Z
Sy=ret.Hy-ret.Hy_y; % Information Storage Y
Txzy=ret.Hy_y-ret.Hy_yzx; % Joint Information Transfer XZ->Y
Txy=ret.Hy_y-ret.Hy_yx; % Information Transfer X->Y
Tzy=ret.Hy_y-ret.Hy_yz; % Information Transfer X->Y
Txy_z=ret.Hy_yz-ret.Hy_yzx; % Conditional Information Transfer X->Y|Z
Tzy_x=ret.Hy_yx-ret.Hy_yzx; % Conditional Information Transfer Z->Y|X
Ixzy=Txzy-(Txy+Tzy); % Interaction transfer XZ->Y
% Partial Information Decomposition
Rxzy=min(Txy,Tzy); % Redundant transfer XZ->Y
Uxy=Txy-Rxzy; % Unique transfer X->Y
Uzy=Tzy-Rxzy; % Unique transfer Z->Y
Sxzy=Txzy-Uxy-Uzy-Rxzy; % Synergistic transfer XZ->Y


%% Estimation on a realization of the simulation
Un = mvnrnd(zeros(1,M),Su,N);
Yn=var_filter(Am,Un); % realization 

% embedding vector
m=p*ones(1,M);
tau=ones(1,M);
VL=its_SetLag(m,tau);

% System Information Y
eHy=its_Elin(Yn(:,jj));
% eHy=out.Hy;

% Information Storage Y
out=its_SElin(Yn,jj,VL);
eSy=out.Sy;
pval_eSy=out.p_Sy;

% Joint Information Transfer XZ->Y
out=its_BTElin(Yn,[ii kk],jj,VL);
eUy_xz=out.Hy_xy;
eTxzy=out.Txy; % joint TE
pval_eTxzy=out.p_Txy;

% Information Transfer X->Y
out=its_BTElin(Yn,ii,jj,VL);
eTxy=out.Txy;
pval_eTxy=out.p_Txy;

% Conditional Information Transfer X->Y|Z
out=its_PTElin(Yn,ii,jj,VL);
eTxy_z=out.Txy_z;
pval_eTxy_z=out.p_Txy_z;

% Derive other quantities
eTzy=eTxzy-eTxy_z; % Information Transfer Z->Y
eTzy_x=eTxzy-eTxy; % Conditional Information Transfer Z->Y|X
eIxzy=eTxzy-(eTxy+eTzy); % Interaction transfer XZ->Y
% Partial Information Decomposition
eRxzy=min(eTxy,eTzy); % Redundant transfer XZ->Y
eUxy=eTxy-eRxzy; % Unique transfer X->Y
eUzy=eTzy-eRxzy; % Unique transfer Z->Y
eSxzy=eTxzy-eUxy-eUzy-eRxzy; % Synergistic transfer XZ->Y





%% display and plot
clc;
disp('Target System information:');
disp(['Hy=' num2str(Hy) ' , eHy=' num2str(eHy)]);
disp(' ');
disp('Information Storage:');
disp(['Sy=' num2str(Sy) ' ; eSy=' num2str(eSy) ', p=' num2str(pval_eSy)]);
disp(' ');
disp('Information Transfer X->Y:');
disp(['Txy=' num2str(Txy) ' ; eTxy=' num2str(eTxy) ', p=' num2str(pval_eTxy)]);
disp(' ');
disp('Conditional Information Transfer X->Y|Z:');
disp(['Txy_z=' num2str(Txy_z) ' ; eTxy_z=' num2str(eTxy_z) ', p=' num2str(pval_eTxy_z)]);
disp(' ');
disp('New information:');
disp(['Ny_xz=' num2str(Uy_xz) ' , eNy_xz=' num2str(eUy_xz)]);
% disp(' ');
% disp('Interaction Information Transfer XZ->Y:');
% disp(['Ixzy=' num2str(Ixzy) ' ; eIxzy=' num2str(eIxzy)]);
% disp(' ');
% disp('Unique Information Transfer X->Y:');
% disp(['Ixzy=' num2str(Uxy) ' ; eIxzy=' num2str(eUxy)]);
% disp(' ');
% disp('Unique Information Transfer Z->Y:');
% disp(['Ixzy=' num2str(Uzy) ' ; eIxzy=' num2str(eUzy)]);
% disp(' ');
% disp('Redundant Information Transfer XZ->Y:');
% disp(['Rxzy=' num2str(Rxzy) ' ; eRxzy=' num2str(eRxzy)]);
% disp(' ');
% disp('Synergistic Information Transfer XZ->Y:');
% disp(['Sxzy=' num2str(Sxzy) ' ; eSxzy=' num2str(eSxzy)]);

figure(1);
subplot(3,3,[1 2]); plot(Yn(:,ii)); zoom xon; title('Xn');
subplot(3,3,3); plot(freq,PSDx); xlim([0 0.5]); title('PSD(Xn)');
subplot(3,3,[4 5]); plot(Yn(:,jj)); zoom xon; title('Yn');
subplot(3,3,6); plot(freq,PSDy); xlim([0 0.5]); title('PSD(Yn)');
subplot(3,3,[7 8]); plot(Yn(:,kk)); zoom xon; title('Zn');
subplot(3,3,9); plot(freq,PSDz); xlim([0 0.5]); title('PSD(Zn)');


