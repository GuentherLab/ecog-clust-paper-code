%% Binning Estimation of the Shannon Entropy
% Computes the shannon entropy of the columns of the matrix A

function e=its_Ebin(Y,c,quantizza,base)

[n,M]=size(Y);

% time series uniform quantization (UQ)
if ~exist('base'), base=exp(1); end
if ~exist('quantizza'), quantizza='y'; end
if quantizza=='n'
    Yq=Y;
else
    for m=1:M
        Yq(:,m)=its_quantization(Y(:,m),c)-1;
    end
end
Q=Yq;

%% compute Shannon entropy
% count identical patterns inside the matrix A (cnt)
% [n,M]=size(A);
% Q=A;
cnt=[];
while ~isempty(Q)
    cmp=[];
    for m=1:M
        cmp=[cmp Q(:,m)==Q(1,m)];
    end
    tmp=(sum(cmp,2)==M);
    cnt=[cnt; sum(tmp)];
    Q(tmp,:)=[];
end

%%% Entropy of A
p=cnt./n;
e=0;
for i=1:length(cnt)
   if p(i)== 0
   else
      e=e-p(i)*log(p(i))/log(base);
   end   
end
