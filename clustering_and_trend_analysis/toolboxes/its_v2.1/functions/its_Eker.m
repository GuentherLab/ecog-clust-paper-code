%% kernel Estimation of the  Entropy 
% Computes the entropy of the M-dimensional variable Y (matrix N*M)
% uses step kernel and maximum (Chebyshev) distance
% Note: uses "SampEn -like" entropy approach (exclude self-matches and compute -log of average probability)

function Hy=its_Eker(Y,r,norma)

if ~exist('norma','var') 
    norma='c'; %default chebishev
end

% Y=[(1:10)' (11:20)'];% Y=(1:10)';
N=size(Y,1);

X=Y';
p=nan*ones(N,1);
for n=1:N
   Xtmp=X; Xtmp(:,n)=[]; % exclude self-matches
   if norma=='c'
       dist = max(abs( repmat(Y(n,:)',1,N-1) - Xtmp ),[],1); % maximum norm
   else % euclidian norm - for dimension 1 is equivalent
       dist=nan*ones(1,size(Xtmp,2));
       for q=1:size(Xtmp,2)
           dist(q)=norm(Y(n,:)'-Xtmp(:,q),2);
       end
   end
    
   D = (dist < r);

   p(n) = sum(D)/(N-1);
end

Hy = - log( mean(p) ); % "SampEn -like" version of entropy estimate 

end









