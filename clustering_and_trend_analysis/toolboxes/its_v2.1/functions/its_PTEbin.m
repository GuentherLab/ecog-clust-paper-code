%% PARTIAL TRANSFER ENTROPY WITH BINNING and NON-UNIFORM EMBEDDING
% Computes the transfer entropy from the column ii to column jj of Y N*M matrix data), conditioned to the remaining columns
% quantifies the Transfer Entropy TEx->y|z using the scheme of Faes EMBC2013+IEEE-TBME2014

%%% inputs are:
% Y: N*M matrix of the M time series each having length N
% ii: index of driver variable
% jj: index of target variable 
% V: vector of candidates (to be pre-determined from uniform embedding or using its_NUEbin.m)
% c: n. of quantization levels for binning
% quantizza: set at 'n' to skip quantization (if you pass data already quantized)

function out=its_PTEbin(Y,V,ii,jj,c,quantizza)

[~,M]=size(Y);

% time series uniform quantization (UQ)
if ~exist('quantizza'), quantizza='y'; end
if quantizza=='n'
    Yq=Y;
else
    for m=1:M
        Yq(:,m)=its_quantization(Y(:,m),c)-1;
    end
end

if isempty(V)
    TE=0;
    Hy=its_Ebin(Yq(:,jj),c,'n');
    out.Txy_z=TE;
    out.Hy_yz=Hy;
    out.Hy_xyz=Hy;
    return
end


B=its_buildvectors(Yq,jj,V); %% form the observation matrix
N=size(B,1);
% set subspaces of lower dimension
A=B; A(:,1)=[];
% XYZ=A;
tmp=V(:,1); iYZ=find(tmp~=ii);
YZ=A(:,iYZ);
yYZ=[B(:,1) YZ];

Hy_xyz=its_CEbin(B);
if isempty(YZ)
    Hy_yz=its_Ebin(Yq(:,jj),c,'n');
else
    Hy_yz=its_CEbin(yYZ);
end

TE=Hy_yz-Hy_xyz;

out.Txy_z=TE;
out.Hy_yz=Hy_yz;
out.Hy_xyz=Hy_xyz;

end










