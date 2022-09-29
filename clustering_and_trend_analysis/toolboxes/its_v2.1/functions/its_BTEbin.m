%% bivariate TRANSFER ENTROPY WITH BINNING and NON-UNIFORM EMBEDDING
% Computes the transfer entropy from a multivariate driver variable (according to ii) to a scalar target variable (column jj of Y)
% V= assigned embedding vector

%%% inputs are:
% Y: data matrix (N*M)
% ii: index of driver variable (can be multivariate)
% jj: index of target variable 
% V: vector of candidates (to be pre-determined from uniform embedding or using CEnu_bin.m)
% c: n. of quantization levels for binning
% quantizza: set at 'n' to skip quantization (if you pass data already quantized)

function out=its_BTEbin(Y,V,ii,jj,c,quantizza)

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
    out.Txy=TE;
    out.Hy_y=Hy;
    out.Hy_xy=Hy;
    return
end



B=its_buildvectors(Yq,jj,V);
% set subspaces of lower dimension
A=B(:,2:end);
tmp=V(:,1);
i_Y = tmp==jj;

% M_Y=A(:,i_Y);
% M_yY=[B(:,1) M_Y];

i_XY=i_Y;
for cni=1:length(ii)
    i_XY= i_XY | tmp==ii(cni);
end
M_y=B(:,1);
M_Y=A(:,i_Y); M_yY=[M_y M_Y];
M_XY=A(:,i_XY); M_yXY=[M_y M_XY];

% Hy_xy=its_CEbin(B);
if isempty(M_Y)
    Hy_y=its_Ebin(Yq,c,'n');
else
    Hy_y=its_CEbin(M_yY);
end

if isempty(M_XY)
    Hy_xy=Hy_y;
else
    Hy_xy=its_CEbin(M_yXY);
end

TE=Hy_y-Hy_xy;

out.Txy=TE;
out.Hy_y=Hy_y;
out.Hy_xy=Hy_xy;


end










