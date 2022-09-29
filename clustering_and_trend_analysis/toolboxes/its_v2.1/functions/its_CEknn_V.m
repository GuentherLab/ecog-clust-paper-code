%% k-nearest neighbor Estimation of the Conditional Entropy
% Extension to the multivariate case (target jj and drivers ii)

function out=its_CEknn_V(data,V,ii,jj,k,metric)
% clear;close all;clc;
% data=[[1 2 3 4 5 6 7 8]' [11 22 33 44 55 66 77 88]'];
% V=[1 1; 1 2; 2 0];
% jj=1; ii=2;
% k=2;

if ~exist('metric','var'), metric='maximum'; end

Hy=its_Eknn(data(:,jj),k,metric);

if isempty(V)
    out.Hy=Hy;
    out.Hy_xy=Hy;
    out.Hy_y=Hy;
    out.Hy_x=Hy;
    return
end

% [Ny,M]=size(Y);

B=its_buildvectors(data,jj,V); %% form the observation matrix
N=size(B,1);

% set subspaces of lower dimension
A=B; A(:,1)=[];
XY=A;
tmp=V(:,1);
iX=find(tmp==ii);
X=A(:,iX);
yX=[B(:,1) X];
iY=find(tmp==jj);
Y=A(:,iY);
yY=[B(:,1) Y];

% neighbor search in space of higher dimension - yXY
atria = nn_prepare(B, metric);
[~, distances] = nn_search(B, atria, (1:N)', k, 0);
dd=distances(:,k);

% range searches in subspaces of lower dimension
atriaXY = nn_prepare(XY, metric); %XY (past of X and Y)
[countXY, neighbors] = range_search(XY, atriaXY, (1:N)', dd, 0);
tmp=neighbors(:,2);%subtraction from count of points with distance exactly equal to k-th neighbor
for n=1:length(tmp)
%     countXYZ(n)=countXYZ(n)-sum(tmp{n}==dd(n));
    countXY(n)=max(k-1,countXY(n)-sum(tmp{n}==dd(n)));
end

atriay = nn_prepare(B(:,1), metric); % y (present of Y)
[county, tmp] = range_search(B(:,1), atriay, (1:N)', dd, 0);
tmp=tmp(:,2); %subtraction from count of points with distance exactly equal to k-th neighbor
for n=1:length(tmp)
%     count_y(n)=count_y(n)-sum(tmp{n}==dd(n));
    county(n)=max(k-1,county(n)-sum(tmp{n}==dd(n)));
end

atriayY = nn_prepare(yY, metric); % y (present of Y)
[countyY, tmp] = range_search(yY, atriayY, (1:N)', dd, 0);
tmp=tmp(:,2); %subtraction from count of points with distance exactly equal to k-th neighbor
for n=1:length(tmp)
%     count_y(n)=count_y(n)-sum(tmp{n}==dd(n));
    countyY(n)=max(k-1,countyY(n)-sum(tmp{n}==dd(n)));
end

atriayX = nn_prepare(yX, metric); % y (present of Y)
[countyX, tmp] = range_search(yX, atriayX, (1:N)', dd, 0);
tmp=tmp(:,2); %subtraction from count of points with distance exactly equal to k-th neighbor
for n=1:length(tmp)
%     count_y(n)=count_y(n)-sum(tmp{n}==dd(n));
    countyX(n)=max(k-1,countyX(n)-sum(tmp{n}==dd(n)));
end


if ~isempty(X) 
    atriaX = nn_prepare(X, metric); % X (past of X)
    [countX, neighbors] = range_search(X, atriaX, (1:N)', dd, 0);
    tmp=neighbors(:,2); %subtraction from count of points with distance exactly equal to k-th neighbor
    for n=1:length(tmp)
%         countYZ(n)=countYZ(n)-sum(tmp{n}==dd(n));
        countX(n)=max(k-1,countX(n)-sum(tmp{n}==dd(n)));
    end
else
    countX=(N-1)*ones(N,1);
end

if ~isempty(Y) 
    atriaY = nn_prepare(Y, metric); % Y (past of Y)
    [countY, neighbors] = range_search(Y, atriaY, (1:N)', dd, 0);
    tmp=neighbors(:,2); %subtraction from count of points with distance exactly equal to k-th neighbor
    for n=1:length(tmp)
%         countYZ(n)=countYZ(n)-sum(tmp{n}==dd(n));
        countY(n)=max(k-1,countY(n)-sum(tmp{n}==dd(n)));
    end
else
    countY=(N-1)*ones(N,1);
end



% equations for conditional entropies
dd2=2*dd;
dd2(dd2==0)=[]; % do not accept distance=0

H_yXY= - psi(k) + psi(N) + (1+length(iX)+length(iY))*mean(log(dd2));
H_XY= psi(N) - (1/N)*(sum(psi(countXY+1))) + (length(iX)+length(iY))*mean(log(dd2));
H_yY= psi(N) - (1/N)*(sum(psi(countyY+1))) + (1+length(iY))*mean(log(dd2));
H_y= psi(N) - (1/N)*(sum(psi(county+1))) + 1*mean(log(dd2));
H_X= psi(N) - (1/N)*(sum(psi(countX+1))) + length(iX)*mean(log(dd2));
H_yX= psi(N) - (1/N)*(sum(psi(countyX+1))) + (1+length(iX))*mean(log(dd2));
H_Y= psi(N) - (1/N)*(sum(psi(countY+1))) + length(iY)*mean(log(dd2));

Hy_XY=H_yXY-H_XY;
Hy_Y=H_yY-H_Y;
Hy_X=H_yX-H_X;

out.Hy=H_y;
out.Hy_xy=Hy_XY;
out.Hy_y=Hy_Y;
out.Hy_x=Hy_X;

out.counts.yYZ=countXY;
out.counts.YZ=countX;
out.counts.XYZ=countY;


    
end

%  figure(11);plot(psi(0.1:0.01:1));

    
