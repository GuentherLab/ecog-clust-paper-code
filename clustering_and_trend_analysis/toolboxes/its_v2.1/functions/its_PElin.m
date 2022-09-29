%% PREDICTIVE INFORMATION WITH LINEAR APPROX and UNIFORM EMBEDDING
% Given a multivariate system of M variables (N*M matrix data), and assigned:
% a target variable Y (j-th col of data)
% a driver variable X (i-th col of data)
% all remaining variables Z (all other columns of data, indexes in vector k)
% quantifies the CE terms useful to compute Eqs. 17(a-d) in Faes,Porta, bookchapter 2013

% inputs:
% Y: N*M data matrix
% V: vector of component indices
% ii: index of driver variable
% jj: index of target variable

function [out]=its_PElin(Y,V,ii,jj)

B=its_buildvectors(Y,jj,V);
A=B(:,2:end);
tmp=V(:,1);
i_Y = tmp==jj;
i_X = tmp==ii;
M=size(Y,2);
kk=setdiff((1:M),[jj ii]); % indexes of z
i_Z = tmp==kk;
i_YZ=i_Y;
i_XZ=i_X;
for cni=1:length(kk)
    i_YZ= i_YZ | tmp==kk(cni);
    i_XZ= i_XZ | tmp==kk(cni);
end
i_XY=i_Y;
for cni=1:length(ii)
    i_XY= i_XY | tmp==ii(cni);
end
M_yXYZ=B;
M_y=B(:,1);
M_XYZ=A;
M_Y=A(:,i_Y); M_yY=[M_y M_Y];
M_X=A(:,i_X); M_yX=[M_y M_X];
M_Z=A(:,i_Z); M_yZ=[M_y M_Z];
M_YZ=A(:,i_YZ); M_yYZ=[M_y M_YZ];
M_XZ=A(:,i_XZ); M_yXZ=[M_y M_XZ];
M_XY=A(:,i_XY); M_yXY=[M_y M_XY];


%% ANALYSIS
[Hy, Ssy]=its_Elin(Y(:,jj));

Uy=Y(:,jj); % residuals of unconditioned regression are the series itself

%%% Unrestricted regression (full-conditioned)
[Hy_xyz,Sy_xyz,Uy_xyz]=its_CElin(M_yXYZ);
% [~,M]=size(Y); zerolag=[1 1 0];
% p_vett=max(V(:,2))*ones(1,M)-zerolag;
% V2=its_SetLag(p_vett,ones(1,M),ones(1,M),zerolag);
% M_yXYZ2=its_buildvectors(Y,jj,V2);
% [Hy_xyz,Sy_xyz,Uy_xyz]=its_CElin(M_yXYZ2);

%%% Restricted regressions (partially-conditioned)

% remove input
if isempty(M_YZ)
    Hy_yz=Hy; Sy_yz=Ssy; Uy_yz=Uy;
else
    [Hy_yz,Sy_yz,Uy_yz]=its_CElin(M_yYZ);
end

% remove target
if isempty(M_XZ)
    Hy_xz=Hy; Sy_xz=Ssy; Uy_xz=Uy;
else
    [Hy_xz,Sy_xz,Uy_xz]=its_CElin(M_yXZ);
end

% remove other
if isempty(M_XY)
    Hy_xy=Hy; Sy_xy=Ssy; Uy_xy=Uy;
else
    [Hy_xy,Sy_xy,Uy_xy]=its_CElin(M_yXY);
end

% use only input
if isempty(M_X)
    Hy_x=Hy; Sy_x=Ssy; Uy_x=Uy;
else
    [Hy_x,Sy_x,Uy_x]=its_CElin(M_yX);
end


% use only target
if isempty(M_Y)
    Hy_y=Hy; Sy_y=Ssy; Uy_y=Uy;
else
    [Hy_y,Sy_y,Uy_y]=its_CElin(M_yY);
end

% use only other
if isempty(M_Z)
    Hy_z=Hy; Sy_z=Ssy; Uy_z=Uy;
else
    [Hy_z,Sy_z,Uy_z]=its_CElin(M_yZ);
end


%%% CondEnt differences and F-tests
% System Predictive Information
Py=Hy-Hy_xyz; 
pval_Py=its_LinReg_Ftest(Uy_xyz,Uy,size(M_XYZ,2),0);
% Self Entropy
Sy=Hy-Hy_y; 
pval_Sy=its_LinReg_Ftest(Uy_y,Uy,size(M_Y,2),0);
% Bivariate Transfer Entropy for z
Tzy=Hy_y-Hy_yz; 
pval_Tzy=its_LinReg_Ftest(Uy_yz,Uy_y,size(M_YZ,2),size(M_Y,2));
% Multivariate Transfer Entropy for x
Txy_z=Hy_yz-Hy_xyz; 
pval_Txy_z=its_LinReg_Ftest(Uy_xyz,Uy_yz,size(M_XYZ,2),size(M_YZ,2));
% Collective Transfer Entropy
Txzy=Hy_y-Hy_xyz; 
pval_Txzy=its_LinReg_Ftest(Uy_xyz,Uy_y,size(M_XYZ,2),size(M_Y,2));

% Bivariate Transfer Entropy for x
Txy=Hy_y-Hy_xy; 
pval_Txy=its_LinReg_Ftest(Uy_xy,Uy_y,size(M_XY,2),size(M_Y,2));
% Multivariate Transfer Entropy for z
Tzy_x=Hy_xy-Hy_xyz; 
pval_Tzy_x=its_LinReg_Ftest(Uy_xyz,Uy_xy,size(M_XYZ,2),size(M_XY,2));
% Interaction Transfer Entropy
Ixzy=Txy-Txy_z;
%Bxz2=Tzy-Tzy_x
%Bxz3=Txy+Tzy-Txzy

%%% arrange output structure
out.Py=Py;
out.Sy=Sy;
out.Tzy=Tzy;
out.Txy_z=Txy_z;
out.Txy=Txy;
out.Tzy_x=Tzy_x;
out.Txzy=Txzy;

out.p_Py=pval_Py;
out.p_Sy=pval_Sy;
out.p_Tzy=pval_Tzy;
out.p_Txy=pval_Txy;
out.p_Txy_z=pval_Txy_z;
out.p_Tzy_x=pval_Tzy_x;
out.p_Txzy=pval_Txzy;

out.H.Hy=Hy;
out.H.Hy_x=Hy_x;
out.H.Hy_y=Hy_y;
out.H.Hy_z=Hy_z;
out.H.Hy_xy=Hy_xy;
out.H.Hy_xz=Hy_xz;
out.H.Hy_yz=Hy_yz;
out.H.Hy_xyz=Hy_xyz;
out.Ixzy=Ixzy;

out.SigmaY=[Ssy Sy_x Sy_y Sy_z Sy_xy Sy_xz Sy_yz Sy_xyz];



%% disp
% disp(['TE' int2str(ii) '->' int2str(jj) '=' num2str(TElin)]);
% disp(['p=' num2str(p_value)]);
% 
% figure(1);
% subplot(3,1,1); plot(Y(:,1)); ylabel('y1');
% subplot(3,1,2); plot(Y(:,2)); ylabel('y2');
% subplot(3,1,3); plot(Y(:,3)); ylabel('y3');

