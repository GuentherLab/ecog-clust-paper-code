%% k-nearest neighbor Estimation of the Self Entropy (Information Storage)
% Computes the information stored in the jj column of Y 
% V= assigned embedding vector (should contain only components of Y, i.e. with jj as first column; anyway others are discarded)

function out=its_SEknn(Y,V,jj,k,metric)

if ~exist('metric','var'), metric='maximum'; end

if isempty(V)
    out.Sy=0;
    out.Hy_y=nan;
    return
end
%% form the observation matrices
B=its_buildvectors(Y,jj,V);
A=B(:,2:end);
tmp=V(:,1);

i_Y= tmp==jj;

M_y=B(:,1);
M_Y=A(:,i_Y);
M_yY=[M_y M_Y];

N=size(B,1);


%% kNN analysis
%%% neighbor search in space of higher dimension
atria_yY = nn_prepare(M_yY, metric);
[~, distances] = nn_search(M_yY, atria_yY, (1:N)', k, 0);
dd=distances(:,k);

%%% range searches in subspaces of lower dimension - M_y
atria_y = nn_prepare(M_y, metric);
[count_y, tmp] = range_search(M_y, atria_y, (1:N)', dd, 0);
tmp=tmp(:,2); %subtraction from count of points with distance exactly equal to k-th neighbor
for n=1:length(tmp)
%     count_y(n)=count_y(n)-sum(tmp{n}==dd(n));
    count_y(n)=max(k-1,count_y(n)-sum(tmp{n}==dd(n)));
end

%%% range searches in subspaces of lower dimension - M_Y
if ~isempty(M_Y)
    atria_Y = nn_prepare(M_Y, metric);
    [count_Y, tmp] = range_search(M_Y, atria_Y, (1:N)', dd, 0);
    tmp=tmp(:,2);%subtraction from count of points with distance exactly equal to k-th neighbor
    for n=1:length(tmp)
%         count_Y(n)=count_Y(n)-sum(tmp{n}==dd(n));
        count_Y(n)=max(k-1,count_Y(n)-sum(tmp{n}==dd(n)));
    end
else
    count_Y=(N-1)*ones(N,1);
end


%% compute PI and terms
Sy = psi(N) + psi(k) -(1/N)*( sum(psi(count_y+1)) + sum(psi(count_Y+1)) );

dd2=2*dd;
dd2(dd2==0)=[]; % do not accept distance=0
Hy_y= - psi(k) + (1/N)*( sum(psi(count_Y+1)) ) + mean(log(dd2)) ;


out.Sy=Sy;
out.Hy_y=Hy_y;

end









