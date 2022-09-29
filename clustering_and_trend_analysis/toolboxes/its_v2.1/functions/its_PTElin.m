%% PARTIAL TANSFER ENTROPY for LINEAR GAUSSIAN PROCESSES (UNIFORM EMBEDDING)

%%% INPUTS:
% data: N*M matrix of the M time series each having length N
% ii: index of driver variable
% jj: index of target variable 
% V: vector of candidates (to be pre-determined from uniform embedding or using its_NUEbin.m)

function out=its_PTElin(data,ii,jj,V)

if isempty(V) % no conditioning: TE=0, CE=E
    [Hy_xyz,Sy_xyz]=its_Elin(data(:,jj));
    out.Hy_yz=Hy_xyz; out.Hy_xyz=Hy_xyz; 
    out.SigmaY_yz=Sy_xyz; out.SigmaY_xyz=Sy_xyz;
    out.Txy_z=0;
    return
end
[~,M]=size(data);

% form the observation matrix
yXYZ=its_buildvectors(data,jj,V);

% set subspaces of lower dimension
XYZ=yXYZ; XYZ(:,1)=[];
tmp=V(:,1); iYZ=find(tmp~=ii);
YZ=XYZ(:,iYZ);
yYZ=[yXYZ(:,1) YZ];

iXZ=find(tmp~=jj);
XZ=XYZ(:,iXZ);
yXZ=[yXYZ(:,1) XZ];

%%% Unrestricted regression (full-conditioned)
[Hy_xyz,Sy_xyz,Uy_xyz,Am]=its_CElin(yXYZ);

%%% Restricted regression (partially-conditioned)
[Hy_yz,Sy_yz,Uy_yz]=its_CElin(yYZ);

TE=Hy_yz-Hy_xyz;
pval = its_LinReg_Ftest(Uy_xyz',Uy_yz',size(V,1),size(YZ,2));

%%% Cross regression
[Hy_xz,Sy_xz]=its_CElin(yXZ);
SE=Hy_xz-Hy_xyz;


%%% output structure
out.Txy_z=TE;
out.p_Txy_z=pval;
out.Hy_yz=Hy_yz;
out.Hy_xyz=Hy_xyz;
out.SigmaY_yz=Sy_yz;
out.SigmaY_xyz=Sy_xyz;
out.Uy_yz=Uy_yz;
out.Uy_xyz=Uy_xyz;
out.coeff=Am;

out.Hy_xz=Hy_xz;
out.SigmaY_xz=Sy_xz;
out.Sy_xz=SE;

end




