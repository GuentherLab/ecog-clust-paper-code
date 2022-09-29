%% LINEAR ANALYSIS OF BRAIN-BODY INTERACTIONS IN REST/STRESS STATES
%% INFORMATION STORAGE AND INFORMATION TRANSFER (total+conditional) for individual nodes at given spatial location and experimental condition
clc; clear; close all;
NumSubj = 18;

% 1   2  3  4   5  6  7  8  9  10 11  12 13 14
% AF3 F7 F3 FC5 T7 P7 O1 O2 P8 T8 FC6 F4 F8 AF4
eegChannel = 1; %spatial location
state=2; % exp.condition: 1)rest 2)mental stress 3)serious game

%%% preprocessing parameters
pfilter=0.94;
p=8; % maximum lag (model order)


%% Load file series and perform network analysis
load('BrainBodyStress.mat');

a=series{1}(:,1:7)
save('stressdata.txt','a','-double')
% fid=fopen('file1.txt','w');
% b=fprintf(fid,'%.5f\n',a)
% fclose(fid)


pepe

for i_subj = 1:NumSubj
    clc; disp(['subject ' int2str(i_subj) '...']);
    data = series{i_subj,state};
    data = data(:,[1:3 4*eegChannel:4*eegChannel+3]);
    [N,M]=size(data);
    
    % filter and normalize to unit variance
    for m=1:M
        dataf(:,m)=AR_filter(data,m,pfilter);
        datan(:,m)=(dataf(:,m)-mean(dataf(:,m)))/std(dataf(:,m));
    end
    Y=datan; %data matrix
        
    % embedding vector
    pV=p*ones(1,M); tau=ones(1,M);
    VL=its_SetLag(pV,tau);
    
    % Storage, total Transfer, New Information
    for jj=1:M
        %%% Information Storage
        out=its_SElin(Y,jj,VL);
        Sy(jj,i_subj)=out.Sy;
        
        %%% Total Information Transfer
        kk=setdiff((1:M),jj);
        out=its_BTElin(Y,kk,jj,VL);
        U{i_subj}(:,jj)=out.Uy_xy'; %(residuals)
        Ny(jj,i_subj)=out.Hy_xy;
        Ty(jj,i_subj)=out.Txy; % joint TE from kk to jj
    end
    
    % conditional information transfer
    for jj=1:M
        for ii=1:M
            if ii~=jj
                out=its_PTElin(Y,ii,jj,VL);    
                Txy(jj,ii,i_subj)=out.Txy_z;
            end
        end
    end
    
end


%% mean values and visualization
Sy_m=mean(Sy,2);
Ty_m=mean(Ty,2);
Ny_m=mean(Ny,2);
Txy_m=mean(Txy,3);

figure(1)
boxplot(Sy');
title('Information Storage');
ylabel('S_Y')
xticklabels({'\eta','\rho','\pi','\delta','\theta','\alhpa','\beta'})

figure(2)
boxplot(Ty');
title('Information Transfer (total)');
ylabel('T_Y')
xticklabels({'\eta','\rho','\pi','\delta','\theta','\alhpa','\beta'})

figure(4)
imagesc(Txy_m); colormap(jet); colorbar
title('Conditional Information Transfer');
xticklabels({'\eta','\rho','\pi','\delta','\theta','\alpha','\beta'})
yticklabels({'\eta','\rho','\pi','\delta','\theta','\alpha','\beta'})

