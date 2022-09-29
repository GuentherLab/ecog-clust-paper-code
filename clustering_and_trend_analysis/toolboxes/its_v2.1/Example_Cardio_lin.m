%% Example of PTE computed with kNN on cardiovascular variability series
clear;close all;clc;

%% Parameters
jj=1; % index of target series (Y: 1=RR)
ii=3; % indexes of input series (X: 2=SAP 3=RESP)

% TElin parameters
pcrit='b'; %model order or criterion 'a' for AIC, 'b' for BIC 'c' for imposed
pmin=1; pmax=20; % max for 'a' and 'b'
pimp=5; % if pcrit='c', choose p=pimp 

%% open cardiovascular data
data=load('supine.prn');
rr=data(:,1);
sap=data(:,2);
resp=data(:,3);
Yo=[rr sap resp];
[N,M]=size(Yo);
kk=setdiff(1:M,[jj ii]); %index of Z

% Remove mean
for m=1:M 
    Y(:,m)=Yo(:,m)-mean(Yo(:,m));
    %Y(:,m)=(Yo(:,m)-mean(Yo(:,m)))/std(Yo(:,m)); 
end


%% Analysis: test optimal model order
zerolag=[0 1 1]; % inst.eff from SAP and resp to RR 
[pottbic,pottaic]=its_FindOrderLin(Yo,jj,pmin,pmax,ones(1,M),ones(1,M),zerolag);
if strcmp(pcrit,'a'), p=pottaic; end
if strcmp(pcrit,'b'), p=pottbic; end
if strcmp(pcrit,'c'), p=pimp; end


%% Information Dynamics
clc;
pV=p*ones(1,M);
VL=its_SetLag(pV,ones(1,M),ones(1,M),zerolag);

% System Information
Hy=its_Elin(Y(:,jj));

% Information Storage
out=its_SElin(Y,jj,VL);
Ny=out.Hy_y;
Sy=out.Sy;

% Information Transfer X->Y
out=its_BTElin(Y,ii,jj,VL);
Txy=out.Txy;
Ny_x=out.Hy_xy;

% Information Transfer XZ->Y
out=its_BTElin(Y,kk,jj,VL);
Ny_xz=out.Hy_xy;
Txzy=out.Txy; % joint TE

% Conditional Information Transfer X->Y|Z
out=its_PTElin(Y,ii,jj,VL);
Hy_yz=out.Hy_yz;
Txy_z=out.Txy_z;

% Information Decomposition
% yields the same results of the above analyses :)
% out=its_PElin(Y,VL,ii,jj);

%% plots and disps
label{1}='RR'; label{2}='SAP'; label{3}='RESP';
figure(1);
subplot(3,1,1); plot(Yo(:,jj)); zoom xon; ylabel(['series ' num2str(jj)]); title(['Target Y: ' label{jj}]);
subplot(3,1,2); plot(Yo(:,ii)); zoom xon; ylabel(['series ' num2str(ii)]); title(['Driver X: ' label{ii}]);
subplot(3,1,3); plot(Yo(:,kk)); zoom xon; ylabel(['series ' num2str(kk)]); title(['Other Z: ' label{kk}]);

clc;
disp(['LINEAR ESTIMATOR']);
disp(['Information Storage ' label{jj} ':']);
disp(['Sy=' num2str(Sy)]); disp(' ');
disp(['Joint Information Transfer from ' label{ii} ',' label{kk} ' to ' label{jj} ':']);
disp(['Txz->y=' num2str(Txzy)]); disp(' ');
disp(['Information Transfer from ' label{ii} ' to ' label{jj} ':']);
disp(['Tx->y=' num2str(Txy)]); disp(' ');
disp(['Conditional Information Transfer from ' label{ii} ' to ' label{jj} ' given ' label{kk} ':']);
disp(['Tx->y|z=' num2str(Txy_z)]); disp(' ');

%%
% % PTE X->Y|Z
% [PTE,p_valPTE,Hy_xyz,Hy_yz,Sy_xyz,Sy_yz,Uy_xyz,Am]=its_PTElin(Yo,ii,jj,VL);
% 
% % GC X->Y|Z (verification - note: this assumes no instantaneous effects!)
% [GC,p_valGC,Am_ver,Su,SigmaR,Ures]=egc_gcMVAR(Yo',p);
% PTE_ver=0.5*GC(jj,ii); p_valPTE_ver=p_valGC(jj,ii);
% 
% %% (b) joint TE
% % JTE XZ->Y
% [JTE,p_valJTE,Hy_xyJ,Hy_yJ,Sy_xyJ,Sy_yJ,Uy_xyJ,AmJ]=its_BTElin(data,jj,VL);
% 
% 
% %% (c) bivariate TE
% kk=setdiff(1:M,[jj ii]); %index of Z
% pV2=p*ones(1,M); pV2(kk)=0;
% zerolag2=zerolag; zerolag2(kk)=0;
% VL2=its_SetLag(pV2,ones(1,M),ones(1,M),zerolag2);
% % BTE X->Y
% [BTE,p_valBTE,Hy_xy,Hy_y,Sy_xy,Sy_y,Uy_xy,Am2]=its_BTElin(data,jj,VL2);
%  
% 
% %% (d) Predictive Information and its decomposition
% out=its_PElin(data,VL,ii,jj);
% % verifications
% PTE_2=out.PI(4); p_valPTE_2=out.pval(4);
% BTE_2=out.PI(5); p_valBTE_2=out.pval(5);
% JTE_2=out.PI(7); p_valJTE_2=out.pval(7);
% 
% %% Plots and disps   
% figure(1); labelfig1='ZYX';
% for m=1:M
%     subplot(M,1,m);
%     plot(Yo(:,m),'.-');
%     xlabel('n'); ylabel([labelfig1(m) '(n)']);
% end
% 
% clc;
% disp(['PTE X->Y|Z=' num2str(PTE)]);
% disp(['JTE XZ->Y=' num2str(JTE)]);
% disp(['BTE X->Y=' num2str(BTE)]);

