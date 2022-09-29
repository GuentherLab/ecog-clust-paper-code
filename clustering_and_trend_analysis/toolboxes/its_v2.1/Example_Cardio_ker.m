%% Computation of Information Dynamics cardiovascular variability series
% KERNEL estimation approach
clear;close all;clc;

% NUembed parameters
rfact=0.3; % kernel threshold set as percent of SD
Lmax=5; %NUembed: number of lags for each system
num_rnd=100; minshift=20;alpha_rnd=0.05; %bin-based NUembed: parameters of termination criterion based on time-shifted surrogates
quantize='y'; %set at 'n' to skip quantization (if you pass data already quantized)

jj=1; % index of target series (Y: 1=RR)
ii=2; % index of input series (X: 2=SAP, or 3=RESP)

%% open cardiovascular data
data=load('supine.prn');
rr=data(:,1);
sap=data(:,2);
resp=data(:,3);
Yo=[rr sap resp];
[N,M]=size(Yo);
kk=setdiff(1:M,[jj ii]); %index of Z

% Normalization - Set 0 mean and 1 stdev
for m=1:M 
    Y(:,m)=(Yo(:,m)-mean(Yo(:,m)))/std(Yo(:,m)); 
end
tau=ones(1,M); %embedding lag always unitary

%% Univariate analysis - system Y
% System Information
r=rfact*std(Y(:,jj));
Hy=its_Eker(Y(:,jj),r);

% Non-uniform Embedding
pV1=Lmax*ones(1,M); pV1(ii)=0; pV1(kk)=0;
zerolag1=zeros(1,M); %set instantaneous effects - here absent (monovariate embedding)
candidates1=its_SetLag(pV1,tau,ones(1,M),zerolag1);
clc; disp('Non uniform embedding - embedding of Y past...');
ret1=its_NUEker(Y,jj,candidates1,r,num_rnd,minshift,alpha_rnd);
VL1=ret1.VL; % embedding vector

% Information Storage, New Information (Univariate)
out=its_SEker(Y,VL1,jj,r);
Sy=out.Sy;
Ny=out.Hy_y;


%% Bivariate analysis - driver:  system X; target: system Y
% Non-uniform Embedding
pV2=Lmax*ones(1,M); pV2(kk)=0;
zerolag2=zeros(1,M); zerolag2(ii)=1; %allow instantaneous effects from X to Y 
candidates2=its_SetLag(pV2,tau,ones(1,M),zerolag2);
clc; disp('Non uniform embedding - embedding of Ypast, Xpast...');
ret2=its_NUEker(Y,jj,candidates2,r,num_rnd,minshift,alpha_rnd);
VL2=ret2.VL; % embedding vector

% Information Transfer, new Information (Bivariate)
out=its_BTEker(Y,VL2,ii,jj,r);
Txy=out.Txy;
Ny_x=out.Hy_xy;


%% Multivariate analysis - drivers:  system X,Z; target: system Y
% Non-uniform Embedding
pV3=Lmax*ones(1,M);
zerolag3=zeros(1,M); zerolag3(ii)=1; zerolag3(kk)=1; %allow instantaneous effects from X to Y and from Z to Y
candidates3=its_SetLag(pV3,tau,ones(1,M),zerolag3);
clc; disp('Non uniform embedding - embedding of Ypast, Xpast, Zpast...');
ret3=its_NUEker(Y,jj,candidates3,r,num_rnd,minshift,alpha_rnd);
VL3=ret3.VL; % embedding vector

% Conditional Information Transfer X->Y|Z, New Information (Multivariate)
out=its_PTEker(Y,VL3,ii,jj,r);
Txy_z=out.Txy_z;
Ny_xz=out.Hy_xyz;

% Joint Information Transfer XZ->Y
out=its_BTEker(Y,VL3,[ii kk],jj,r);
Txzy=out.Txy;


%% plots and disps
label{1}='RR'; label{2}='SAP'; label{3}='RESP';
figure(1);
subplot(3,1,1); plot(Yo(:,jj)); zoom xon; ylabel(['series ' num2str(jj)]); title(['Target Y: ' label{jj}]);
subplot(3,1,2); plot(Yo(:,ii)); zoom xon; ylabel(['series ' num2str(ii)]); title(['Driver X: ' label{ii}]);
subplot(3,1,3); plot(Yo(:,kk)); zoom xon; ylabel(['series ' num2str(kk)]); title(['Other Z: ' label{kk}]);

clc;
disp(['KERNEL ESTIMATOR']);
disp(['Information Storage ' label{jj} ':']);
disp(['Sy=' num2str(Sy)]); disp(' ');
disp(['Joint Information Transfer from ' label{ii} ',' label{kk} ' to ' label{jj} ':']);
disp(['Txz->y=' num2str(Txzy)]); disp(' ');
disp(['Information Transfer from ' label{ii} ' to ' label{jj} ':']);
disp(['Tx->y=' num2str(Txy)]); disp(' ');
disp(['Conditional Information Transfer from ' label{ii} ' to ' label{jj} ' given ' label{kk} ':']);
disp(['Tx->y|z=' num2str(Txy_z)]); disp(' ');




