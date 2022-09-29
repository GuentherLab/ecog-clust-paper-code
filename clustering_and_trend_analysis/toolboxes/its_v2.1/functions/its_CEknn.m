%% k-nearest neighbor Estimation of the Conditional Entropy
% Computes the conditional information of the first column of B conditioned to the remaining columns (A)
% k is the number of neighbors, 'metric' works for maximum (Chebyshev) distance for now (Euclidian distance should also be implemented

function Hy_y=its_CEknn(B,k,metric)

if ~exist('metric','var'), metric='maximum'; end

A=B(:,2:end);
N=size(B,1);


%% kNN analysis
%%% neighbor search in space of higher dimension
atria_yY = nn_prepare(B, metric);
[~, distances] = nn_search(B, atria_yY, (1:N)', k, 0);
dd=distances(:,k);

%%% range search in the subspace of lower dimension - M_Y
if ~isempty(A)
    atria_Y = nn_prepare(A, metric);
    [count_Y, tmp] = range_search(A, atria_Y, (1:N)', dd, 0);
    tmp=tmp(:,2);%subtraction from count of points with distance exactly equal to k-th neighbor
    for n=1:length(tmp)
%         count_Y(n)=count_Y(n)-sum(tmp{n}==dd(n));
        count_Y(n)=max(k-1,count_Y(n)-sum(tmp{n}==dd(n)));
    end
else
    count_Y=(N-1)*ones(N,1);
end


%% computes CE

dd2=2*dd;
dd2(dd2==0)=[]; % do not accept distance=0
Hy_y= - psi(k) + (1/N)*( sum(psi(count_Y+1)) ) + mean(log(dd2)) ;


end









