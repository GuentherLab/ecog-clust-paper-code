%% k-nearest neighbor Estimation of the Bivariate Transfer Entropy
% Computes the bivariate transfer entropy from a multivariate driver variable (according to ii) to a scalar target variable (column jj of Y)
% V= assigned embedding vector

function out=its_BTEknn(Y,V,ii,jj,k,metric)

if ~exist('metric','var'), metric='maximum'; end

if isempty(V)
    out.Txy=0;
    out.Hy_xy=nan;
    return
end
%% form the observation matrices
B=its_buildvectors(Y,jj,V);
A=B(:,2:end);
tmp=V(:,1);

i_Y= tmp==jj;
M_Y=A(:,i_Y);

% M_yY=[B(:,1) M_Y];
% M_YZ=A;
% M_yYZ=B;

i_YZ=i_Y;
for cni=1:length(ii)
    i_YZ= i_YZ | tmp==ii(cni);
end
M_y=B(:,1);
M_yY=[M_y M_Y];
M_YZ=A(:,i_YZ);
M_yYZ=[M_y M_YZ];

%[Ny,M]=size(Y);
N=size(B,1);


%% kNN analysis
%%% neighbor search in space of higher dimension
atria_yYZ = nn_prepare(M_yYZ, metric);
[~, distances] = nn_search(M_yYZ, atria_yYZ, (1:N)', k, 0);
dd=distances(:,k);

%%% range searches in subspaces of lower dimension - M_yY
atria_yY = nn_prepare(M_yY, metric);
[count_yY, tmp] = range_search(M_yY, atria_yY, (1:N)', dd, 0);
tmp=tmp(:,2); %subtraction from count of points with distance exactly equal to k-th neighbor
for n=1:length(tmp)
%     count_yY(n)=count_yY(n)-sum(tmp{n}==dd(n));
    count_yY(n)=max(k-1,count_yY(n)-sum(tmp{n}==dd(n)));
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

%%% range searches in subspaces of lower dimension - M_YZ
if ~isempty(M_YZ)
    atria_YZ = nn_prepare(M_YZ, metric);
    [count_YZ, tmp] = range_search(M_YZ, atria_YZ, (1:N)', dd, 0);
    tmp=tmp(:,2);%subtraction from count of points with distance exactly equal to k-th neighbor
    for n=1:length(tmp)
%         count_YZ(n)=count_YZ(n)-sum(tmp{n}==dd(n));
        count_YZ(n)=max(k-1,count_YZ(n)-sum(tmp{n}==dd(n)));
    end
else
    count_YZ=(N-1)*ones(N,1);
end


%% compute BTE
TE = psi(k) + (1/N)*( sum(psi(count_Y+1)) - sum(psi(count_yY+1)) - sum(psi(count_YZ+1)) );

dd2=2*dd;
dd2(dd2==0)=[]; % do not accept distance=0
Hy_yz= - psi(k) + (1/N)*( sum(psi(count_YZ+1))) + mean(log(dd2));

out.Txy=TE;
out.Hy_xy=Hy_yz;

end
