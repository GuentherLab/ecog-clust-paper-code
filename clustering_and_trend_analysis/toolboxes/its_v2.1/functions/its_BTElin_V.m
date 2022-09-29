%% BIVARIATE TRANSFER ENTROPY for LINEAR GAUSSIAN PROCESSES (UNIFORM EMBEDDING)
%% version for vector target process: computes TE from ii to jj
%%% INPUTS:
% data: N*M matrix of the M time series each having length N
% ii: index of driver variable (can be multivariate)
% jj: index of target variable (can be multivariate)
% p: vector of candidates (to be pre-determined from uniform embedding)

function out=its_BTElin_V(data,ii,jj,p)

X=data(:,ii); [N,P]=size(X);
Y=data(:,jj); Q=size(Y,2);

Z=[Y X];
Vr=its_SetLag(p*ones(1,Q),ones(1,Q));
Vu=its_SetLag(p*ones(1,Q+P),ones(1,Q+P));

for j=1:Q
    B=its_buildvectors(Z,j,Vr);
    [~,~,Uj]=its_CElin(B);
    Ur(:,j)=Uj'; 
    
    B=its_buildvectors(Z,j,Vu);
    [~,~,Uj]=its_CElin(B);
    Uu(:,j)=Uj';
    
end
SigmaUr=cov(Ur);
SigmaUu=cov(Uu);

HUr=0.5*log(det(SigmaUr))+0.5*Q*log(2*pi*exp(1));
HUu=0.5*log(det(SigmaUu))+0.5*Q*log(2*pi*exp(1));


% ce=0.5*log(det(S))+0.5*size(Yb,1)*log(2*pi*exp(1)); %e.g., Barnett PRL 2009
TE=HUr-HUu;

pval = its_LinReg_Ftest_var(det(SigmaUu),det(SigmaUr),p*(P+Q),p*Q,N);

out.TE=TE;
out.pval_TE=pval;
out.SigmaUr=SigmaUr;
out.SigmaUu=SigmaUu;

end
