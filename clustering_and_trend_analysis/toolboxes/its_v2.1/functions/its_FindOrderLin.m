%% Finds the optimal model order for univariate regressions using AIC or BIC reasonings
% works on a single series to be predicted from itself and other M-1 series
% considers also possible instantaneous predictions (through 'zerolag')

% inputs:
% Y: data, N observations (rows) and M variables (columns)
% jj: index of target series (to be predicted)
% pmin,, pmax: range of searches for optimal order
% tau: vector of embedding delays (one for each series)
% u: vector of propagation times (one for each series)
% zerolag: for each series, 1 if zerolag effect is wanted, 0 if not

function [p_bic,p_aic,aic,bic]=its_FindOrderLin(Y,jj,pmin,pmax,tau,u,zerolag)

[N,M]=size(Y);

Var_jj=NaN*ones(pmax,1);
aic=NaN*ones(pmax,1); bic=NaN*ones(pmax,1);
for ip=pmin:pmax
    p_tmp=ip*ones(1,M)-zerolag; % for each series, ip terms in total
    V_tmp=its_SetLag(p_tmp,tau,u,zerolag);
    B=its_buildvectors(Y,jj,V_tmp);
    if isempty(V_tmp) % no conditioning: CE=E
        [~,Var_jj(ip)]=its_Elin(Y(:,jj));
    end
    [~,Var_jj(ip)]=its_CElin(B);
    k=size(V_tmp,1);
    aic(ip)=N*log(Var_jj(ip))+2*k;
    bic(ip)=N*log(Var_jj(ip))+log(N)*k;
end

p_aic=find(aic == min(aic));
p_bic=find(bic == min(bic));

