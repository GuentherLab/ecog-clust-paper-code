%% k-nearest neighbor Estimation of the Mutual Information
%% between two Vector processes
% Computes the mutual information between two vector variables (identified with indexes ii and jj)
% Z= data matrix (series in column)
% ii, jj: row vectors containing the indexes of the columns of Z corresponding to the two variables X and Y for which to compute MI
% k: numer of neighbors for MI computation

%%%%%%%%%%%%%%%TEST
% clear;close all;clc;
% Z=randn(10,7); %10 data points, 7 series
% ii=[5 6 7]; jj=[1 2];
% k=2; metric='maximum';

function MI=its_MIknn_V(Z,ii,jj,k,metric)
if ~exist('metric','var'), metric='maximum'; end

%% form the observation matrices
dimjj=length(jj);
dimii=length(ii);
V=[ii',zeros(dimii,1)];
M_YX=its_buildvectors(Z,jj,V);
M_Y=M_YX(:,1:dimjj);
M_X=M_YX(:,dimjj+1:end);

N=size(M_YX,1);


%% kNN analysis
%%% neighbor search in space of higher dimension - M_YX
atria_YX = nn_prepare(M_YX, metric);
[~, distances] = nn_search(M_YX, atria_YX, (1:N)', k, 0);
dd=distances(:,k);

%%% range searches in subspaces of lower dimension - M_Y
atria_Y = nn_prepare(M_Y, metric);
[count_Y, tmp] = range_search(M_Y, atria_Y, (1:N)', dd, 0);
tmp=tmp(:,2); %subtraction from count of points with distance exactly equal to k-th neighbor
for n=1:length(tmp)
%     count_yY(n)=count_yY(n)-sum(tmp{n}==dd(n));
    count_Y(n)=max(k-1,count_Y(n)-sum(tmp{n}==dd(n)));
end

%%% range searches in subspaces of lower dimension - M_X
atria_X = nn_prepare(M_X, metric);
[count_X, tmp] = range_search(M_X, atria_X, (1:N)', dd, 0);
tmp=tmp(:,2); %subtraction from count of points with distance exactly equal to k-th neighbor
for n=1:length(tmp)
%     count_yY(n)=count_yY(n)-sum(tmp{n}==dd(n));
    count_X(n)=max(k-1,count_X(n)-sum(tmp{n}==dd(n)));
end



%% compute MI
MI = psi(N) + psi(k) - (1/N)*( sum(psi(count_X+1))+sum(psi(count_Y+1)) ) ;

end
