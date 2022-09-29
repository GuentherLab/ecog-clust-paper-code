%% computation of Conditional entropy from VAR model parameters
% Version for monovariate time series (computes Entropy and Storage)

% Instead of using the data, uses the VAR parameters (either estimated from data, or theoretical)
% Works under the linear Gaussian assumption


function ret = its_CElinVAR1(Am,Su,q)


%% PART 1: COVARIANCE MATRICES
M = size(Am,1); %number of elements in system
p=floor(size(Am,2)/M); %number of lags in MVAR model

R=NaN*ones(M,M,q+1); % prepare covariance matrices, (:,:,1) is lag 0, (:,:,q+1) is lag q 

% Obtain F and Delta
Im=eye(M*p);
F=[Am;Im(1:end-M,:)];% this is A^p
Delta=zeros(p*M,p*M); %(this is actually Sigma^p, but use Delta for clarity in code)
Delta(1:M,1:M)=Su(:,:);

% Obtain R_o^p=BigSigma solving the Lyapunov equation: BigSigma = F * BigSigma * F^T + Delta
BigSigma=dlyap(F,Delta);

% extract R(0),...,R(p-1)
for i=1:p
    R(:,:,i)=BigSigma(1:M,M*(i-1)+1:M*i);
end

% Yule-Walker solution  for lags >= p
for k=p+1:q+1
    Rk=R(:,:,k-1:-1:k-p);
    Rm=[];
    for ki=1:p
        Rm=[Rm; Rk(:,:,ki)];
    end
    R(:,:,k)=Am*Rm;
end




%% PART 2: CONDITIONAL ENTROPIES
L=q; % maximum lag of Xpast, Ypast, Zpast

Ry=NaN*ones(1,q+1);
for k=1:q+1
    Ry(k)=R(1,1,k);
end

%%% Entropy of Y
Sy=Ry(1);
Hy=0.5*log(Sy)+0.5*log(2*pi*exp(1));

%%% Conditional Entropy of Y given Ypast
SYpast=toeplitz(Ry(1:L),Ry(1:L));
SyYpast=Ry(2:L+1);
% Sy_Ypast = Sy - SyYpast * inv(SYpast) * SyYpast';
Sy_Ypast = Sy - SyYpast/SYpast * SyYpast';
Hy_Ypast=0.5*log(Sy_Ypast)+0.5*log(2*pi*exp(1));



ret.Sy=Sy;
ret.Sy_y=Sy_Ypast;
ret.Hy=Hy;
ret.Hy_y=Hy_Ypast;

ret.BigSigma=BigSigma;
ret.R=R;

end

