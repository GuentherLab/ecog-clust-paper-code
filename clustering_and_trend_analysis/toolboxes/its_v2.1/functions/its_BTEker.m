%% bivariate TRANSFER ENTROPY WITH KERNEL and NON-UNIFORM EMBEDDING
% Computes the transfer entropy from a multivariate driver variable (according to ii) to a scalar target variable (column jj of Y)
% V= assigned embedding vector

%%% inputs are:
% Y: data matrix (N*M)
% ii: index of driver variable (can be multivariate)
% jj: index of target variable 
% V: vector of candidates (to be pre-determined from uniform embedding or using CEnu_bin.m)
% r is threshold for distances, pass the absolute value
% norma: type of distance ('c' for Chebyshev distance, euclidean otherwise)

function out=its_BTEker(Y,V,ii,jj,r,norma)

if ~exist('norma','var') 
    norma='c'; %default Chebyshev
end

[~,M]=size(Y);

if isempty(V)
    TE=0;
    Hy=its_Eker(Y(:,jj),r,norma);
    out.Txy=TE;
    out.Hy_y=Hy;
    out.Hy_xy=Hy;
    return
end



B=its_buildvectors(Y,jj,V);
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
    Hy_y=its_Eker(Y(:,jj),r,norma);
else
    Hy_y=its_CEker(M_yY,r,norma);
end

if isempty(M_XY)
    Hy_xy=Hy_y;
else
    Hy_xy=its_CEker(M_yXY,r,norma);
end

TE=Hy_y-Hy_xy;

out.Txy=TE;
out.Hy_y=Hy_y;
out.Hy_xy=Hy_xy;


end










