%% PARTIAL TRANSFER ENTROPY WITH KERNEL and NON-UNIFORM EMBEDDING
% Computes the transfer entropy from the column ii to column jj of Y N*M matrix data), conditioned to the remaining columns
% quantifies the Transfer Entropy TEx->y|z using the scheme of Faes EMBC2013+IEEE-TBME2014

%%% inputs are:
% Y: N*M matrix of the M time series each having length N
% ii: index of driver variable
% jj: index of target variable 
% V: vector of candidates (to be pre-determined from uniform embedding or using its_NUEbin.m)
% r is threshold for distances, pass the absolute value
% norma: type of distance ('c' for Chebyshev distance, euclidean otherwise)

function out=its_PTEker(Y,V,ii,jj,r,norma)

if ~exist('norma','var') 
    norma='c'; %default Chebyshev
end

%[~,M]=size(Y);

if isempty(V)
    TE=0;
    Hy=its_Eker(Y(:,jj),r,norma);
    out.Txy_z=TE;
    out.Hy_yz=Hy;
    out.Hy_xyz=Hy;
    return
end


B=its_buildvectors(Y,jj,V); %% form the observation matrix
%N=size(B,1);
% set subspaces of lower dimension
A=B; A(:,1)=[];
% XYZ=A;
tmp=V(:,1); iYZ=find(tmp~=ii);
YZ=A(:,iYZ);
yYZ=[B(:,1) YZ];

Hy_xyz=its_CEker(B,r,norma);
if isempty(YZ)
    Hy_yz=its_Eker(Y(:,jj),r,norma);
else
    Hy_yz=its_CEker(yYZ,r,norma);
end

TE=Hy_yz-Hy_xyz;

out.Txy_z=TE;
out.Hy_yz=Hy_yz;
out.Hy_xyz=Hy_xyz;

end










