%% k-nearest neighbor Estimation of the Bivariate Transfer Entropy from two sources to a target + individual TEs
%% from Vector source process to Scalar target process
% Computes the bivariate transfer entropy from two multivariate driver variables (according to ii and kk)
% to a scalar target variable (column jj of Y)
% Y= data matrix (series in column)
% V= assigned embedding vector

function out=its_JBTEknn_VS(Y,V,ii,kk,jj,knn,metric)

if ~exist('metric','var'), metric='maximum'; end

% V=[1 1; 1 2; 2 1; 3 1; 3 2; 4 1; 5 1];
% Y=rand(10,5);
% jj=3; %i_Y
% ii=[1 4]; %i_X1
% kk=2; % i_X2


if isempty(V)
    out.Tx1y=0;
    out.Tx2y=0;
    out.Tx1x2y=0;
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

iikk=[ii kk];
i_YX1X2=i_Y;
for cni=1:length(iikk)
    i_YX1X2= i_YX1X2 | tmp==iikk(cni);
end
i_YX1=i_Y;
for cni=1:length(ii)
    i_YX1= i_YX1 | tmp==ii(cni);
end
i_YX2=i_Y;
for cni=1:length(kk)
    i_YX2= i_YX2 | tmp==kk(cni);
end


M_y=B(:,1);
M_yY=[M_y M_Y];
M_YX1X2=A(:,i_YX1X2);
M_yYX1X2=[M_y M_YX1X2];
M_YX1=A(:,i_YX1);
M_yYX1=[M_y M_YX1];
M_YX2=A(:,i_YX2);
M_yYX2=[M_y M_YX2];

%[Ny,M]=size(Y);
N=size(B,1);


%% kNN analysis
%%% neighbor search in space of higher dimension - M_yYX1X2
atria_yYZ = nn_prepare(M_yYX1X2, metric);
[~, distances] = nn_search(M_yYX1X2, atria_yYZ, (1:N)', knn, 0);
dd=distances(:,knn);

%%% range searches in subspaces of lower dimension - M_yY
atria_yY = nn_prepare(M_yY, metric);
[count_yY, tmp] = range_search(M_yY, atria_yY, (1:N)', dd, 0);
tmp=tmp(:,2); %subtraction from count of points with distance exactly equal to k-th neighbor
for n=1:length(tmp)
%     count_yY(n)=count_yY(n)-sum(tmp{n}==dd(n));
    count_yY(n)=max(knn-1,count_yY(n)-sum(tmp{n}==dd(n)));
end

%%% range searches in subspaces of lower dimension - M_Y
if ~isempty(M_Y)
    atria_Y = nn_prepare(M_Y, metric);
    [count_Y, tmp] = range_search(M_Y, atria_Y, (1:N)', dd, 0);
    tmp=tmp(:,2);%subtraction from count of points with distance exactly equal to k-th neighbor
    for n=1:length(tmp)
%         count_Y(n)=count_Y(n)-sum(tmp{n}==dd(n));
        count_Y(n)=max(knn-1,count_Y(n)-sum(tmp{n}==dd(n)));
    end
else
    count_Y=(N-1)*ones(N,1);
end

%%% range searches in subspaces of lower dimension - M_YX1X2
if ~isempty(M_YX1X2)
    atria_YX1X2 = nn_prepare(M_YX1X2, metric);
    [count_YX1X2, tmp] = range_search(M_YX1X2, atria_YX1X2, (1:N)', dd, 0);
    tmp=tmp(:,2);%subtraction from count of points with distance exactly equal to k-th neighbor
    for n=1:length(tmp)
%         count_YZ(n)=count_YZ(n)-sum(tmp{n}==dd(n));
        count_YX1X2(n)=max(knn-1,count_YX1X2(n)-sum(tmp{n}==dd(n)));
    end
else
    count_YX1X2=(N-1)*ones(N,1);
end

%%% range searches in subspaces of lower dimension - M_YX1
if ~isempty(M_YX1)
    atria_YX1 = nn_prepare(M_YX1, metric);
    [count_YX1, tmp] = range_search(M_YX1, atria_YX1, (1:N)', dd, 0);
    tmp=tmp(:,2);%subtraction from count of points with distance exactly equal to k-th neighbor
    for n=1:length(tmp)
%         count_YZ(n)=count_YZ(n)-sum(tmp{n}==dd(n));
        count_YX1(n)=max(knn-1,count_YX1(n)-sum(tmp{n}==dd(n)));
    end
else
    count_YX1=(N-1)*ones(N,1);
end

%%% range searches in subspaces of lower dimension - M_yYX1
atria_yYX1 = nn_prepare(M_yYX1, metric);
[count_yYX1, tmp] = range_search(M_yYX1, atria_yYX1, (1:N)', dd, 0);
tmp=tmp(:,2); %subtraction from count of points with distance exactly equal to k-th neighbor
for n=1:length(tmp)
%     count_yY(n)=count_yY(n)-sum(tmp{n}==dd(n));
    count_yYX1(n)=max(knn-1,count_yYX1(n)-sum(tmp{n}==dd(n)));
end

%%% range searches in subspaces of lower dimension - M_YX2
if ~isempty(M_YX2)
    atria_YX2 = nn_prepare(M_YX2, metric);
    [count_YX2, tmp] = range_search(M_YX2, atria_YX2, (1:N)', dd, 0);
    tmp=tmp(:,2);%subtraction from count of points with distance exactly equal to k-th neighbor
    for n=1:length(tmp)
%         count_YZ(n)=count_YZ(n)-sum(tmp{n}==dd(n));
        count_YX2(n)=max(knn-1,count_YX2(n)-sum(tmp{n}==dd(n)));
    end
else
    count_YX2=(N-1)*ones(N,1);
end

%%% range searches in subspaces of lower dimension - M_yYX2
atria_yYX2 = nn_prepare(M_yYX2, metric);
[count_yYX2, tmp] = range_search(M_yYX2, atria_yYX2, (1:N)', dd, 0);
tmp=tmp(:,2); %subtraction from count of points with distance exactly equal to k-th neighbor
for n=1:length(tmp)
%     count_yY(n)=count_yY(n)-sum(tmp{n}==dd(n));
    count_yYX2(n)=max(knn-1,count_yYX2(n)-sum(tmp{n}==dd(n)));
end



%% compute joint TE X1X2->Y and individual TEs X1->y and X2->Y
TEX1X2_Y = psi(knn) + (1/N)*( sum(psi(count_Y+1)) - sum(psi(count_yY+1)) - sum(psi(count_YX1X2+1)) );

TEX1_Y = (1/N)*( sum(psi(count_Y+1)) + sum(psi(count_yYX1+1)) - sum(psi(count_yY+1)) - sum(psi(count_YX1+1)) );

TEX2_Y = (1/N)*( sum(psi(count_Y+1)) + sum(psi(count_yYX2+1)) - sum(psi(count_yY+1)) - sum(psi(count_YX2+1)) );

TEX1_Y_X2 = psi(knn) + (1/N)*( -sum(psi(count_yYX2+1)) + sum(psi(count_YX2+1)) - sum(psi(count_YX1X2+1)) );

TEX2_Y_X1 = psi(knn) + (1/N)*( -sum(psi(count_yYX1+1)) + sum(psi(count_YX1+1)) - sum(psi(count_YX1X2+1)) );


dd2=2*dd;
dd2(dd2==0)=[]; % do not accept distance=0
Hy_yx1x2= - psi(knn) + (1/N)*( sum(psi(count_YX1X2+1))) + mean(log(dd2));

out.Tx1x2y=TEX1X2_Y;
out.Tx1y=TEX1_Y; % TE from ii to jj
out.Tx2y=TEX2_Y; % TE from kk to jj

out.TEX1_Y_X2=TEX1_Y_X2; % TE from ii to jj given kk
out.TEX2_Y_X1=TEX2_Y_X1; % TE from kk to jj given ii
    
out.Hy_x1x2y=Hy_yx1x2;

end
