%% Linear Estimation of the Mutual Information between two Vector processes
% Computes the mutual information between two vector variables (identified with indexes ii and jj)
% Z= data matrix (series in column)

function [MI,pval]=its_MIlin_V(Z,ii,jj)

X=Z(:,ii);
Y=Z(:,jj);
[N,P]=size(X);
Q=size(Y,2);

SigmaX=cov(X);
HX=0.5*log(det(SigmaX))+0.5*P*log(2*pi*exp(1));
for i=1:P
    Ai=[X(:,i) Y];
    [~,~,Ui_oth]=its_CElin(Ai);
    U(:,i)=Ui_oth';
end
SigmaU=cov(U);
HU=0.5*log(det(SigmaU))+0.5*P*log(2*pi*exp(1));
MI=HX-HU;
pval = its_LinReg_Ftest_var(det(SigmaU),det(SigmaX),Q,0,N);



end
