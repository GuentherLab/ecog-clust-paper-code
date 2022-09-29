function ret = its_CElinVAR0(Am,Su,jj,ii,q,zerolag)
% extension of CElin_analytic that accounts for possible zerolag effects

%   Computes Information Dynamics analytically for a stationary mvar(p) process:
%   X_n=A_1*X_{n-1}+A_2*X_{t-n}+...+A_p*X_{n-p}+E_n
%
%   INPUT: 
%   Am  -   generalized connectivity matrix A=(A_1 A_2 ... A_p)
%   Su  -   covariance matrix for E_n
%   jj  -   index of target series (y)
%   ii  -   index of input series (x)
%   q   -   number of lags
%   zerolag- 1x2 vector of flags of instantaneous effects [x z]=[ii kk], with 1 where one wants the zerolag effect
%   OUTPUT:
%   ret structure with all possible conditional entropies


%%%%% internal test
% variables to be passed are Am, Su, q, ii, jj, zerolag
% clear; close all; clc;
% zerolag=[1 1];
% % [Am,Ak,Su,lambdamax,~]=coeffMVAR(0);
% [Am,Ak,Su,lambdamax,Yo]=coeffMVAR4(0);
% q=5; % maximum lag desired for the covariance matrices
% jj=3; % index of target variable Y (scalar)
% ii=1; % index of input variable X (scalar)

%% PART 1: COVARIANCE MATRICES
M = size(Am,1); %number of elements in system
p=floor(size(Am,2)/M); %number of lags in MVAR model

R=NaN*ones(M,M,q+1); % prepare covariance matrices, (:,:,1) is lag 0, (:,:,q+1) is lag q 

% Obtain F and Delta
Im=eye(M*p);
F=[Am;Im(1:end-M,:)];% this is A^p
Delta=zeros(p*M,p*M); %(this is actually Sigma^p, but use Delta for clarity in code)
Delta(1:M,1:M)=Su(:,:);

% Obtain R_o^p=BigSigma solving the Lyapunov equation: BigSigma = F * BigSigma * F^T + Delta
BigSigma=dlyap(F,Delta);

% extract R(0),...,R(p-1)
for i=1:p
    R(:,:,i)=BigSigma(1:M,M*(i-1)+1:M*i);
end

% Yule-Walker solution  for lags >= p
for k=p+1:q+1
    Rk=R(:,:,k-1:-1:k-p);
    Rm=[];
    for ki=1:p
        Rm=[Rm; Rk(:,:,ki)];
    end
    R(:,:,k)=Am*Rm;
end




%% PART 2: CONDITIONAL ENTROPIES
kk=setdiff([1:M],[ii jj]); % index of Z (vector of dimension M-2)
L=q; % maximum lag of Xpast, Ypast, Zpast

%%% create covariance matrices related to X,Y,Z
Ry=NaN*ones(1,q+1);
Rx=NaN*ones(1,q+1);
Rz=NaN*ones(M-2,M-2,q+1);
Rxy=NaN*ones(1,q+1);
Ryx=NaN*ones(1,q+1);
Ryz=NaN*ones(1,M-2,q+1);
Rzy=NaN*ones(M-2,1,q+1);
Rxz=NaN*ones(1,M-2,q+1);
Rzx=NaN*ones(M-2,1,q+1);
for k=1:q+1
    Ry(k)=R(jj,jj,k);
    Rx(k)=R(ii,ii,k);
    Rz(:,:,k)=R(kk,kk,k);
    Rxy(k)=R(ii,jj,k);
    Ryx(k)=R(jj,ii,k);
    Ryz(:,:,k)=R(jj,kk,k);
    Rzy(:,:,k)=R(kk,jj,k);
    Rxz(:,:,k)=R(ii,kk,k);
    Rzx(:,:,k)=R(kk,ii,k);
end

%% Entropy of Y
Sy=Ry(1);
Hy=0.5*log(Sy)+0.5*log(2*pi*exp(1));


%% Conditional Entropy of Y given Ypast
SYpast=toeplitz(Ry(1:L),Ry(1:L));
SyYpast=Ry(2:L+1);
% Sy_Ypast = Sy - SyYpast * inv(SYpast) * SyYpast';
Sy_Ypast = Sy - SyYpast/SYpast * SyYpast';
Hy_Ypast=0.5*log(Sy_Ypast)+0.5*log(2*pi*exp(1));



%% Conditional Entropy of Y given Ypast and Zpast
SYZ=NaN*ones(L,(M-2)*L);
for ri=1:L
    for co=1:L
        if co>=ri
            SYZ(ri,(co-1)*(M-2)+1:co*(M-2))=Ryz(:,:,co-ri+1);
        else
            SYZ(ri,(co-1)*(M-2)+1:co*(M-2))=Rzy(:,:,ri-co+1)';
        end
    end
end
SZY=NaN*ones((M-2)*L,L);
for ri=1:L
    for co=1:L
        if co>=ri
            SZY((ri-1)*(M-2)+1:ri*(M-2),co)=Rzy(:,:,co-ri+1);
        else
            SZY((ri-1)*(M-2)+1:ri*(M-2),co)=Ryz(:,:,ri-co+1)';
        end
    end
end
SZ=NaN*ones((M-2)*L,(M-2)*L);
for ri=1:L
    for co=1:L
        if co>=ri
            SZ((ri-1)*(M-2)+1:ri*(M-2),(co-1)*(M-2)+1:co*(M-2))=Rz(:,:,co-ri+1);
        else
            SZ((ri-1)*(M-2)+1:ri*(M-2),(co-1)*(M-2)+1:co*(M-2))=Rz(:,:,ri-co+1)';
        end
    end
end

SYpastZpast=[[SYpast SYZ];[SZY SZ]];
tmp=squeeze(Ryz(:,:,2:L+1));
SyZpast=reshape(tmp,1,L*(M-2));
SyYpastZpast=[SyYpast SyZpast];

%%%% eventual further computations if zerolag effects are present from Z to Y
if zerolag(2)==1 % if Z has instantaneous influence on Y
    SZ0Ypast=NaN*ones(M-2,L);
    SYpastZ0=NaN*ones(L,M-2);
    SZ0Zpast=NaN*ones(M-2,(M-2)*L);
    SZpastZ0=NaN*ones((M-2)*L,M-2);
    for co=1:L
        SZ0Ypast(:,co)=Rzy(:,:,co+1);
        SYpastZ0(co,:)=Rzy(:,:,co+1)';
        SZ0Zpast(:,(co-1)*(M-2)+1:co*(M-2))=Rz(:,:,co+1);
        SZpastZ0((co-1)*(M-2)+1:co*(M-2),:)=Rz(:,:,co+1)';
    end
    SYpastZpast=[ [SYpastZpast [SYpastZ0; SZpastZ0]]; [SZ0Ypast SZ0Zpast Rz(:,:,1)] ];
    SyYpastZpast=[SyYpastZpast Ryz(:,:,1)];
    
    %these will be useful later when conditioning also on Xpast
    SZ0Xpast=NaN*ones(M-2,L);
    SXpastZ0=NaN*ones(L,M-2);
    for co=1:L
        SZ0Xpast(:,co)=Rzx(:,:,co+1);
        SXpastZ0(co,:)=Rzx(:,:,co+1)';
    end
end
% Sy_YpastZpast = Sy - SyYpastZpast * inv(SYpastZpast) * SyYpastZpast'; 
Sy_YpastZpast = Sy - SyYpastZpast/SYpastZpast * SyYpastZpast';
Hy_YpastZpast=0.5*log(Sy_YpastZpast)+0.5*log(2*pi*exp(1));



%% Conditional Entropy of Y given Zpast (added 30/9/2014)
if zerolag(2)==1 % if Z has instantaneous influence on Y
    SZ0=[[SZ SZpastZ0]; [SZ0Zpast Rz(:,:,1)]];
    SyZpast0=[SyZpast Ryz(:,:,1)];
else
    SZ0=SZ;
    SyZpast0=SyZpast;
end
Sy_Zpast = Sy - SyZpast0/SZ0 * SyZpast0';
Hy_Zpast=0.5*log(Sy_Zpast)+0.5*log(2*pi*exp(1));


%% Conditional Entropy of Y given Ypast and Xpast (added 4/1/2017)
SX=toeplitz(Rx(1:L),Rx(1:L));
SyXpast=Ryx(2:L+1);
SYX=toeplitz(Rxy(1:L),Ryx(1:L));
SXY=SYX';
SYpastXpast=[[SYpast SYX];[SXY SX]];
SyYpastXpast=[SyYpast SyXpast];

%%%% eventual further computations if zerolag effects are present from Z to Y
if zerolag(1)==1 % if Z has instantaneous influence on Y
    SYpastX0=Rxy(2:L+1)';
    SX0Ypast=SYpastX0';
    SXpastX0=Rx(2:L+1)';
    SX0Xpast=SXpastX0';
    SYpastXpast=[[SYpastXpast [SYpastX0; SXpastX0]] ; [SX0Ypast SX0Xpast Rx(1)]];
    SyYpastXpast=[SyYpastXpast Ryx(1)];
end
Sy_YpastXpast = Sy - SyYpastXpast/SYpastXpast * SyYpastXpast';
Hy_YpastXpast=0.5*log(Sy_YpastXpast)+0.5*log(2*pi*exp(1));





%% Conditional Entropy of Y given Ypast, Zpast and Xpast
SX=toeplitz(Rx(1:L),Rx(1:L));
SYX=toeplitz(Rxy(1:L),Ryx(1:L));
SXY=SYX';
SXZ=NaN*ones(L,(M-2)*L);
for ri=1:L
    for co=1:L
        if co>=ri
            SXZ(ri,(co-1)*(M-2)+1:co*(M-2))=Rxz(:,:,co-ri+1);
        else
            SXZ(ri,(co-1)*(M-2)+1:co*(M-2))=Rzx(:,:,ri-co+1)';
        end
    end
end
SZX=NaN*ones((M-2)*L,L);
for ri=1:L
    for co=1:L
        if co>=ri
            SZX((ri-1)*(M-2)+1:ri*(M-2),co)=Rzx(:,:,co-ri+1);
        else
            SZX((ri-1)*(M-2)+1:ri*(M-2),co)=Rxz(:,:,ri-co+1)';
        end
    end
end
SYpastZpastXpast=[[SYpast SYZ SYX];[SZY SZ SZX];[SXY SXZ SX]];
SyXpast=Ryx(2:L+1);
SyYpastZpastXpast=[SyYpast SyZpast SyXpast];

%%%% eventual matrix extensions if zerolag effects are present
if zerolag(1)==1 % if X has instantaneous influence on Y
    SX0Ypast=Rxy(2:L+1);
    SYpastX0=Rxy(2:L+1)';
    SX0Xpast=Rx(2:L+1);
    SXpastX0=Rx(2:L+1)';
    SX0Zpast=NaN*ones(1,(M-2)*L);
    SZpastX0=NaN*ones((M-2)*L,1);
    for co=1:L
        SX0Zpast(:,(co-1)*(M-2)+1:co*(M-2))=Rxz(:,:,co+1);
        SZpastX0((co-1)*(M-2)+1:co*(M-2),:)=Rxz(:,:,co+1)';
    end
end

if zerolag(2)==1 % if Z has instantaneous influence on Y
    if zerolag(1)==1 % if also X has instantaneous influence on Y
        SYpastZpastXpast=[ [SYpastZpastXpast [SYpastZ0; SZpastZ0; SXpastZ0] [SYpastX0; SZpastX0; SXpastX0]]; [[SZ0Ypast SZ0Zpast SZ0Xpast Rz(:,:,1) Rzx(:,:,1)];[SX0Ypast SX0Zpast SX0Xpast Rxz(:,:,1) Rx(1)]] ];
        SyYpastZpastXpast=[ SyYpastZpastXpast Ryz(:,:,1) Ryx(1)];
    else % Z has instantaneous influence on Y but X has not
        SYpastZpastXpast=[ [SYpastZpastXpast [SYpastZ0; SZpastZ0; SXpastZ0]] ; [SZ0Ypast SZ0Zpast SZ0Xpast Rz(:,:,1)] ];
        SyYpastZpastXpast=[ SyYpastZpastXpast Ryz(:,:,1)];
    end
else
    if zerolag(1)==1 % if X has instantaneous influence on Y, but Z has not
        SYpastZpastXpast=[ [SYpastZpastXpast [SYpastX0; SZpastX0; SXpastX0]] ; [SX0Ypast SX0Zpast SX0Xpast Rx(1)] ];
        SyYpastZpastXpast=[ SyYpastZpastXpast Ryx(1)];
    end
end
% Sy_YpastZpastXpast = Sy - SyYpastZpastXpast * inv(SYpastZpastXpast) * SyYpastZpastXpast';
Sy_YpastZpastXpast = Sy - SyYpastZpastXpast/SYpastZpastXpast * SyYpastZpastXpast';
Hy_YpastZpastXpast=0.5*log(Sy_YpastZpastXpast)+0.5*log(2*pi*exp(1));



%% Conditional Entropy of Y given Zpast and Xpast (added 30/9/2014)
SZpastXpast=[[SZ SZX];[SXZ SX]];
SyZpastXpast=[SyZpast SyXpast];

if zerolag(2)==1 % if Z has instantaneous influence on Y
    if zerolag(1)==1 % if also X has instantaneous influence on Y
        SZpastXpast=[ [SZpastXpast [SZpastZ0; SXpastZ0] [SZpastX0; SXpastX0]]; [SZ0Zpast SZ0Xpast Rz(:,:,1) Rzx(:,:,1)]; [SX0Zpast SX0Xpast Rxz(:,:,1) Rx(1)]  ];
        SyZpastXpast=[SyZpastXpast Ryz(:,:,1) Ryx(1)];
    else % Z has instantaneous influence on Y but X has not
        SZpastXpast=[ [SZpastXpast [SZpastZ0; SXpastZ0]] ; [SZ0Zpast SZ0Xpast Rz(:,:,1)] ];
        SyZpastXpast=[SyZpastXpast Ryz(:,:,1)];
    end
else
    if zerolag(1)==1 % if X has instantaneous influence on Y, but Z has not
        SZpastXpast=[ [SZpastXpast [SZpastX0; SXpastX0]] ; [SX0Zpast SX0Xpast Rx(1)] ];
        SyZpastXpast=[SyZpastXpast Ryx(1)];
    end
end
Sy_ZpastXpast = Sy - SyZpastXpast/SZpastXpast * SyZpastXpast';
Hy_ZpastXpast=0.5*log(Sy_ZpastXpast)+0.5*log(2*pi*exp(1));




%% return
ret.Hy=Hy; ret.Sy=Sy;
ret.Hy_y=Hy_Ypast; ret.Sy_y=Sy_Ypast;
ret.Hy_yz=Hy_YpastZpast; ret.Sy_yz=Sy_YpastZpast; 
ret.Hy_yx=Hy_YpastXpast; ret.Sy_yx=Sy_YpastXpast; 
ret.Hy_yzx=Hy_YpastZpastXpast; ret.Sy_yzx=Sy_YpastZpastXpast;

ret.Hy_z=Hy_Zpast; ret.Sy_z=Sy_Zpast;
ret.Hy_zx=Hy_ZpastXpast; ret.Sy_zx=Sy_ZpastXpast;

ret.BigSigma=BigSigma;

