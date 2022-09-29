%% Application of framework of Information Dynamics to EEG data
clear; close all; clc;

%%%% PARAMETERS
jj=3; % index of target electrode

set_tau='y'; % 'y' to set time delay at decorrelation time
soglia=1/exp(1); maxlags=10;
filtra='y'; pfilter=0.94; % 'y' to apply detrending filter

% NUembed parameters
knn=10; % n. of neighbors
metric='maximum'; % metric (use always maximum !!)
Lmax=5; %NUembed: number of lags for each system

% file to open
nome='EEG_EyesClosed.mat';

%% load EEG data
load([nome]);

data=v_salvato.ser;
can=v_salvato.nomec';
fc=v_salvato.fc;

data=data(:,12:30);
can=can(12:30); % EEG channels (10-20 system)

[N,M]=size(data);

%%%% re-referencing: Common Average Reference
Y=nan*ones(N,M);
for m=1:M
    tmp=data; tmp(:,m)=[];
    avgY=mean(tmp')';
    Y(:,m)=data(:,m)-avgY;
end

%%%% possible detreding (high-pass filter for each channel)
if filtra=='y'
    for m=1:M
        Y(:,m)=AR_filter(Y,m,pfilter);
    end
end

%%% evaluation of embedding delays (for each series)
if set_tau=='y'
    for i=1:M
        Z=Y(:,i);
        [zc,lags]=xcov(Z,maxlags);
        zcc=zc./zc(maxlags+1); % normalized autocorrelation
        zCC(:,i)=zcc(maxlags+1:2*maxlags+1);
        k=1;
        while k<maxlags+1
            if zCC(k,i)>soglia && zCC(k+1,i)<=soglia
                ztau=k;
                break;
            else
                k=k+1;
            end
        end
        if k==maxlags+1
            tau(i)=maxlags;
        else
            tau(i)=ztau;
        end
    end
else
    tau=ones(1,M); % take all delays identically equal to 1
end


%% Information domain analysis for target jj (system Y)

%%% System Information %%%
Hy=its_Eknn(Y(:,jj),knn);


%%% STORAGE %%%
% Candidates: Lmax lagged components, spaced in time by tau(jj)
pV1=zeros(1,M); pV1(jj)=Lmax;
zerolag1=zeros(1,M); u=ones(1,M); u(jj)=tau(jj);
candidates1=its_SetLag(pV1,tau,u,zerolag1);

% Non-uniform Embedding
clc; disp(['Channel ' int2str(jj) ': ' char(can(jj))]);
disp(' '); disp('NUE Y past...');
minshift=0; numshift=100; alphashift=0.05; % CMI-based NUembed: surro-based termination criterion 
ret1=its_NUEknn(Y,jj,candidates1,knn,numshift,minshift,alphashift);
VL1=ret1.VL; % embedding vector
out=its_SEknn(Y,VL1,jj,knn,metric);

% Information Storage
Sy=out.Sy; % SE



%%% TRANSFER %%%
% Candidates: Lmax lagged components for each channel ch, spaced in time by tau(ch)
pV2=Lmax*ones(1,M); pV2(jj)=0;
zerolag2=zeros(1,M);
u2=ones(1,M);
candidates2=its_SetLag(pV2,tau,u2,zerolag2);
% Non-uniform Embedding   
disp(' ');
disp('NUE X past Z past...');
minshift=max(candidates2(:,2)); numshift=100; alphashift=0.05;
ret2=its_NUEknn_Vstart(Y,jj,candidates2,knn,numshift,minshift,alphashift,metric,VL1); %incremental NUE (starts from V1)
VL2=ret2.VL; % embedding vector

% Global information transfer X->Y, where X is all drivers
ii=VL2(:,1); ii(ii==jj)=[]; ii=unique(ii);
out=its_BTEknn(Y,VL2,ii,jj,knn,metric);
TXy=out.Txy; % TE

% Conditional Information Transfer Xi->Y|Z
Txy_z=zeros(M,1);
for iis=1:length(ii)
    out=its_PTEknn(Y,VL2,ii(iis),jj,k,metric);
    Txy_z(ii(iis))=out.Txy_z; % PTE
end


%% Test for instantaneous effects
insteff=nan*ones(M,1);
disp(' ');
for iik=1:M
    disp(['insteff, input ' num2str(iik)]);
    if iik~=jj
    B=its_buildvectors(Y,jj,[VL2; [iik 0]]);
    CMI=its_CMIknn(B,knn,metric);
    maxshift=size(B,1)-minshift;
    for is=1:numshift
        lagshift=fix(rand*(maxshift-minshift+1)+minshift);% random shift in [minshift, maxshift] 
        xs=circshift(B(:,size(B,2)),lagshift);%is-th shift of the values of the selcted candidate component
        lagshift1=fix(rand*(maxshift-minshift+1)+minshift);
        x1s=circshift(B(:,1),lagshift1);
        Bs=B; Bs(:,end)=xs; % replace last column of B with the shifted term
        Bs(:,1)=x1s; % randomize also target variable
        CMIs(is)=its_CMIknn(Bs,knn,metric);
    end
    soglia=prctile(CMIs,100*(1-alphashift));
    insteff(iik)=CMI>soglia;
    end
end

%% Information domain analysis for target jj (system Y) - COMPENSATED ANALYSIS

% self-embedding vector enlarged with instantaneous effects
tmp=find(insteff==1);
VL1_0=[VL1; [tmp zeros(size(tmp,1),1)]];

%%% TRANSFER %%%
disp(' ');
disp('Compensated TE Analysis...');
disp('NUE X past Z past...');
minshift=0; % minshift=max(candidates2(:,2));
ret3=its_NUEknn_Vstart(Y,jj,candidates2,knn,numshift,minshift,alphashift,metric,VL1_0);
VL3=ret3.VL; % embedding vector

% Global COMPENSATED information transfer X->Y, where X is all drivers
VL3tmp=VL3; VL3tmp(VL3tmp(:,2)==0,:)=[];
ii=VL3tmp(:,1); ii(ii==jj)=[]; ii=unique(ii);
out=its_BTEknn_0(Y,VL3,ii,jj,knn,metric);
cTXy=out.Txy; % cTE

% Conditional COMPENSATED Information Transfer Xi->Y|Z
cTxy_z=zeros(M,1);
for iis=1:length(ii)
    out=its_PTEknn_0(Y,VL3,ii(iis),jj,k,metric);
    cTxy_z(ii(iis))=out.Txy_z; % cPTE
end

%% display results
clc;
disp(['Information Dynamics for channel ' char(can(jj))]);
disp(' ');
disp('Information Storage:');
disp(['SE=' num2str(Sy)]);
disp(' ');
disp('Joint Information Transfer (original, compensated):');
disp(['TE=' num2str(TXy) ', cTE=' num2str(cTXy)]);
disp(' ');
disp('Conditional Information Transfer');
disp(['Channels with significant PTE directed to ' char(can(jj)) ':']);
disp(can(Txy_z>0));
disp(['Channels with significant cPTE directed to ' char(can(jj)) ':']);
disp(can(cTxy_z>0));






