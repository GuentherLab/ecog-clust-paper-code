%% form observation matrix (for entropy computation)
% Y: (quantized) input multiple time series, dimension N*M
% V: list of candidates, dimension Nc*2, Nc is number of candidates; 1st column: index of the signal; 2nd column: index of the lag
% A: output matrix of the vectors specified from the signals Y according to the list V
% B: complete matrix with added the current samples as first column

function B=its_buildvectors(Y,j,V)

% clear;close all;clc;
% tmp=[5 0 0 0 3 5 4 3 4 2 1 2 4 0 5 5]'; Y =[tmp tmp tmp];
% Y=[(1:10)' (11:20)' (21:30)'];
% V=[1 1; 1 3; 2 4; 3 1; 2 1]; % V=[y1(n-1),y1(n-3),y2(n-4),y3(n-1),y2(n-1)]
% V=[1 1; 2 2];
% j=1;


if isempty(V) %if no conditioning, give back the jth series of Y
    B=Y(:,j);
else
    [N,M]=size(Y);
    Nc=size(V,1); % number of candidates

    Lmax=max(V(:,2)); %maximum lag (across all signals)

    A=NaN*ones(N-Lmax,Nc);
    for n=Lmax+1:N
        for i=1:Nc %fill the i-th row of A
            A(n-Lmax,i)=Y(n-V(i,2),V(i,1)); 
        end
    end
    B=[Y((Lmax+1:N)',j) A]; % add current value
end
