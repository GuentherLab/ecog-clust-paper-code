%% Linear Gaussian Estimation of the Shannon Entropy
% Computes the shannon entropy of a multivariate dataset A
% A is N*M multivariate data (N observations, M variables)

function [e, covA]=its_Elin(A)


%% compute Conditional entropy
covA=cov(A);

%Entropy for the multivariate Gaussian case:
e=0.5*log(det(covA))+0.5*size(A,2)*log(2*pi*exp(1)); %e.g., Barnett PRL 2009


% out.Hy=e;
% out.SigmaY=covA;

end

% 0.5*log(det(cov(Q)))+0.5*log(2*pi*exp(1)); 
% log(std(Q)*sqrt(2*pi*exp(1)))