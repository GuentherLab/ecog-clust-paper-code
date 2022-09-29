%% SELF ENTROPY WITH KERNEL ESTIMATOR
% Computes the self entropy of a scalar target variable (column jj of Y)
% V= assigned embedding vector (only the components of jj are used)
% this uses maximum distance (Chebyshev) norm (Euclidian distance should also be implemented)

%%% inputs are:
% Y: data matrix (N*M)
% jj: index of target variable 
% V: vector of candidates (to be pre-determined from uniform embedding or using CEnu_bin.m)
% r is threshold for distances, pass the absolute value
% norma: type of distance ('c' for Chebyshev distance, euclidean otherwise)

function out=its_SEker(Y,V,jj,r,norma)

Y=Y(:,jj); %we will use only the target

if ~exist('norma','var') 
    norma='c'; %default Chebyshev
end

Hy=its_Eker(Y,r,norma);

if isempty(V)
    out.Hy=Hy;
    out.Hy_y=Hy;
    out.Sy=0;
%     SE=0;
    return
end


tmp=V(:,1);
V=V(tmp==jj,:); %remove components other than jj

M_yY=its_buildvectors(Y,1,V);
M_Y=M_yY(:,2:end);

if isempty(M_Y)
    Hy_y=Hy;
else
    Hy_y=its_CEker(M_yY,r,norma);
end

Sy=Hy-Hy_y;

out.Hy=Hy;
out.Hy_y=Hy_y;
out.Sy=Sy;


end










