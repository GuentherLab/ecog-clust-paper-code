%% k-nearest neighbor Estimation of the Predictive Information and decomposition
% Computes Py=Sy+Tzy+Txy_z from an embedding matrix B
% V= assigned embedding vector
% ii= indexes of X; jj=index of Y

function [Py,Sy,Tzy,Txy_z,out]=its_PEknn(Y,V,ii,jj,k,metric)

if ~exist('metric','var'), metric='maximum'; end

if isempty(V)
    Py=0; Sy=0; Tzy=0; Txy_z=0;
    out=NaN;
    return
end
%% form the observation matrices
B=its_buildvectors(Y,jj,V);
A=B(:,2:end);
tmp=V(:,1);

i_Y= tmp==jj;

% i_YZ= tmp~=ii;
% new 12/12/2013: to allow vector-valued ii
M=size(Y,2); kk=setdiff((1:M),[jj ii]);
i_YZ=i_Y;
for cni=1:length(kk)
    i_YZ= i_YZ | tmp==kk(cni);
end


M_yXYZ=B;
M_y=B(:,1);
M_XYZ=A;
M_Y=A(:,i_Y);
M_YZ=A(:,i_YZ);
M_yY=[M_y M_Y];
M_yYZ=[M_y M_YZ];

%[Ny,M]=size(Y);
N=size(B,1);


%% kNN analysis
%%% neighbor search in space of higher dimension
atria_yXYZ = nn_prepare(M_yXYZ, metric);
[~, distances] = nn_search(M_yXYZ, atria_yXYZ, (1:N)', k, 0);
dd=distances(:,k);

%%% range searches in subspaces of lower dimension - M_y
atria_y = nn_prepare(M_y, metric);
[count_y, tmp] = range_search(M_y, atria_y, (1:N)', dd, 0);
tmp=tmp(:,2); %subtraction from count of points with distance exactly equal to k-th neighbor
for n=1:length(tmp)
%     count_y(n)=count_y(n)-sum(tmp{n}==dd(n));
    count_y(n)=max(count_y(n)-sum(tmp{n}==dd(n)),k-1);
end

%%% range searches in subspaces of lower dimension - M_yY
atria_yY = nn_prepare(M_yY, metric);
[count_yY, tmp] = range_search(M_yY, atria_yY, (1:N)', dd, 0);
tmp=tmp(:,2); %subtraction from count of points with distance exactly equal to k-th neighbor
for n=1:length(tmp)
%     count_yY(n)=count_yY(n)-sum(tmp{n}==dd(n));
    count_yY(n)=max(count_yY(n)-sum(tmp{n}==dd(n)),k-1);
end

%%% range searches in subspaces of lower dimension - M_yYZ
atria_yYZ = nn_prepare(M_yYZ, metric);
[count_yYZ, tmp] = range_search(M_yYZ, atria_yYZ, (1:N)', dd, 0);
tmp=tmp(:,2); %subtraction from count of points with distance exactly equal to k-th neighbor
for n=1:length(tmp)
%     count_yYZ(n)=count_yYZ(n)-sum(tmp{n}==dd(n));
    count_yYZ(n)=max(count_yYZ(n)-sum(tmp{n}==dd(n)),k-1);
end

%%% range searches in subspaces of lower dimension - M_Y
if ~isempty(M_Y)
    atria_Y = nn_prepare(M_Y, metric);
    [count_Y, tmp] = range_search(M_Y, atria_Y, (1:N)', dd, 0);
    tmp=tmp(:,2);%subtraction from count of points with distance exactly equal to k-th neighbor
    for n=1:length(tmp)
%         count_Y(n)=count_Y(n)-sum(tmp{n}==dd(n));
        count_Y(n)=max(count_Y(n)-sum(tmp{n}==dd(n)),k-1);
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
        count_YZ(n)=max(count_YZ(n)-sum(tmp{n}==dd(n)),k-1);
    end
else
    count_YZ=(N-1)*ones(N,1);
end

%%% range searche in subspaces of lower dimension - M_XYZ
if ~isempty(M_XYZ) % questo è pleonastico, dato l'if isempty(V) all'inizio...
    atria_XYZ = nn_prepare(M_XYZ, metric);
    [count_XYZ, tmp] = range_search(M_XYZ, atria_XYZ, (1:N)', dd, 0);
    tmp=tmp(:,2);%subtraction from count of points with distance exactly equal to k-th neighbor
    for n=1:length(tmp)
%         count_XYZ(n)=count_XYZ(n)-sum(tmp{n}==dd(n));
        count_XYZ(n)=max(count_XYZ(n)-sum(tmp{n}==dd(n)),k-1);
    end
else
    count_XYZ=(N-1)*ones(N,1);
end




%% compute PI and terms

Py = psi(k) + psi(N) - (1/N)*( sum(psi(count_y+1)) + sum(psi(count_XYZ+1)) );

Sy = psi(N) + (1/N)*( sum(psi(count_yY+1)) - sum(psi(count_y+1)) - sum(psi(count_Y+1)) );

Tzy = (1/N)*( sum(psi(count_Y+1)) + sum(psi(count_yYZ+1)) - sum(psi(count_yY+1)) - sum(psi(count_YZ+1)) );

Txy_z = psi(k) + (1/N)*( sum(psi(count_YZ+1)) - sum(psi(count_yYZ+1)) - sum(psi(count_XYZ+1)) );


out.Ny=count_y;
out.NXYZ=count_XYZ;
out.NyY=count_yY;
out.NY=count_Y;
out.NyYZ=count_yYZ;
out.NYZ=count_YZ;
out.NyYZ=count_yYZ;
out.NXYZ=count_XYZ;







