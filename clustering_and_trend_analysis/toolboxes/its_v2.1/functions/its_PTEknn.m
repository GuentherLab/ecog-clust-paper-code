%% k-nearest neighbor Estimation of the Partial Transfer Entropy
% Computes the transfer entropy from the column ii to column jj of Y, conditioned to the remaining columns

function out=its_PTEknn(Y,V,ii,jj,k,metric)

if ~exist('metric','var'), metric='maximum'; end

if isempty(V)
    TE=0;
    out.Txy_z=TE;
    return
end

% [Ny,M]=size(Y);

B=its_buildvectors(Y,jj,V); %% form the observation matrix
N=size(B,1);

% set subspaces of lower dimension
A=B; A(:,1)=[];
XYZ=A;
tmp=V(:,1); iYZ=find(tmp~=ii);
YZ=A(:,iYZ);
yYZ=[B(:,1) YZ];

% neighbor search in space of higher dimension
atria = nn_prepare(B, metric);
[~, distances] = nn_search(B, atria, (1:N)', k, 0);
dd=distances(:,k);


% range searches in subspaces of lower dimension
atriaXYZ = nn_prepare(XYZ, metric);
[countXYZ, neighbors] = range_search(XYZ, atriaXYZ, (1:N)', dd, 0);
tmp=neighbors(:,2);%subtraction from count of points with distance exactly equal to k-th neighbor
for n=1:length(tmp)
%     countXYZ(n)=countXYZ(n)-sum(tmp{n}==dd(n));
    countXYZ(n)=max(k-1,countXYZ(n)-sum(tmp{n}==dd(n)));
end

atriayYZ = nn_prepare(yYZ, metric);
[countyYZ, neighbors] = range_search(yYZ, atriayYZ, (1:N)', dd, 0);
tmp=neighbors(:,2); %subtraction from count of points with distance exactly equal to k-th neighbor
for n=1:length(tmp)
%     countyYZ(n)=countyYZ(n)-sum(tmp{n}==dd(n));
    countyYZ(n)=max(k-1,countyYZ(n)-sum(tmp{n}==dd(n)));
end


if ~isempty(YZ)
    atriaYZ = nn_prepare(YZ, metric);
    [countYZ, neighbors] = range_search(YZ, atriaYZ, (1:N)', dd, 0);
    tmp=neighbors(:,2); %subtraction from count of points with distance exactly equal to k-th neighbor
    for n=1:length(tmp)
%         countYZ(n)=countYZ(n)-sum(tmp{n}==dd(n));
        countYZ(n)=max(k-1,countYZ(n)-sum(tmp{n}==dd(n)));
    end
else
    countYZ=(N-1)*ones(N,1);
end

% equation of wibral for transfer entropy
TE = psi(k) + (1/N) * ( sum(psi(countYZ+1)) - sum(psi(countyYZ+1)) - sum(psi(countXYZ+1)) ); 

dd2=2*dd;
dd2(dd2==0)=[]; % do not accept distance=0
Hy_xyz= - psi(k) + (1/N)*(sum(psi(countXYZ+1))) + mean(log(dd2));

out.Txy_z=TE;
out.Hy_xyz=Hy_xyz;

out.counts.yYZ=countyYZ;
out.counts.YZ=countYZ;
out.counts.XYZ=countXYZ;
    
end

%  figure(11);plot(psi(0.1:0.01:1));

    
