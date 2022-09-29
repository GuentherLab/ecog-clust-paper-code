%% Linear Gaussian Estimation of the Conditional Entropy
% Computes the conditional entropy of the first column of B conditioned to the remaining columns (A)

% returns estimates of conditional entropy (ce) and residual covariance (S)
% and of the prediction errors (Up, residuals) and the linear regression coefficients (Am)

function [ce,S,Up,Am]=its_CElin(B)


%% compute Conditional entropy
Yb=B(:,1)'; % inversion works with data organized in rows
A=B; A(:,1)=[]; Z=A';

Am=Yb/Z; % least squares!

Yp=Am*Z; 
Up=Yb-Yp;
S=cov(Up');

%CondEn for the Gaussian case:
ce=0.5*log(det(S))+0.5*size(Yb,1)*log(2*pi*exp(1)); %e.g., Barnett PRL 2009
end