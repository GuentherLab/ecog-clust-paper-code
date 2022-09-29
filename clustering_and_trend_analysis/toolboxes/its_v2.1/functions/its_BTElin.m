%% BIVARIATE TRANSFER ENTROPY for LINEAR GAUSSIAN PROCESSES (UNIFORM EMBEDDING)
%%% INPUTS:
% data: N*M matrix of the M time series each having length N
% ii: index of driver variable (can be multivariate)
% jj: index of target variable 
% V: vector of candidates (to be pre-determined from uniform embedding)

function out=its_BTElin(data,ii,jj,V)

if isempty(V) % no conditioning: TE=0, CE=E
    [Hy_xy,Sy_xy]=its_Elin(data(:,jj));
    out.Hy_y=Hy_xy; 
    out.Hy_xy=Hy_xy; 
    out.SigmaY_y=Sy_xy;
    out.SigmaY_xy=Sy_xy;    
    out.Txy=0;
    return
end
% [~,M]=size(data);

% form the observation matrix and set subspaces of lower dimension
B=its_buildvectors(data,jj,V);
A=B(:,2:end);
tmp=V(:,1);
i_Y = tmp==jj;
% i_X = tmp==ii;

i_XY=i_Y; i_X=zeros(size(tmp,1),1);
for cni=1:length(ii)
    i_X=i_X | tmp==ii(cni);
    i_XY= i_XY | tmp==ii(cni);
end
M_y=B(:,1);
M_Y=A(:,i_Y); M_yY=[M_y M_Y];
M_XY=A(:,i_XY); M_yXY=[M_y M_XY];
M_X=A(:,i_X); M_yX=[M_y M_X];

%%% Unrestricted regression (full-conditioned)
[Hy_xy,Sy_xy,Uy_xy,Am]=its_CElin(M_yXY);

%%% Restricted regression (partially-conditioned)
[Hy_y,Sy_y,Uy_y]=its_CElin(M_yY);

TE=Hy_y-Hy_xy;
pval = its_LinReg_Ftest(Uy_xy',Uy_y',size(M_XY,2),size(M_Y,2)); %size(V,1)

%%% Cross regression
[Hy_x,Sy_x]=its_CElin(M_yX);
SE=Hy_x-Hy_xy;

%%% output structure
out.Txy=TE;
out.p_Txy=pval;
out.Hy_y=Hy_y;
out.Hy_xy=Hy_xy;
out.coeff=Am;
out.SigmaY_y=Sy_y;
out.SigmaY_xy=Sy_xy;
out.Uy_xy=Uy_xy;

out.Hy_x=Hy_x;
out.SigmaY_x=Sy_x;
out.Sy_x=SE;

end




