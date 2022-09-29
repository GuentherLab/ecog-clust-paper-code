%% k-nearest neighbor Estimation of the  Entropy 
% Computes the entropy of the M-dimensional variable Y (matrix N*M)

function [Hy,dd]=its_Eknn(Y,k,metric)

if ~exist('metric','var'), metric='maximum'; end

[N,M]=size(Y);

%%% neighbor search in space of higher dimension
atria = nn_prepare(Y, metric);
[~, distances] = nn_search(Y, atria, (1:N)', k, 0);
dd=distances(:,k);

dd(dd==0)=[]; % do not accept distance=0

Hy= - psi(k) + psi(N) + M*mean(log(2*dd));


end









