%% Binning Estimation of the Conditional Entropy
% Computes the conditional entropy of the first column of B conditioned to the remaining columns (A)
% uses the fast version which recalls countmember.m
% NB: B has to be in quantized form

function ce=its_CEbin(B)

%% compute Conditional entropy
yj=B(:,1); 
A=B; A(:,1)=[];

% count identical patterns inside restricted matrix A (cnt) and total matrix B (cnt2)
[n,M]=size(A);
Q=A;
cnt=[];cnt2=[];
while ~isempty(Q)
    cmp=[];
    for m=1:M
        cmp=[cmp Q(:,m)==Q(1,m)];
    end
    tmp=(sum(cmp,2)==M);
    cnt=[cnt; sum(tmp)];
    cnt2=[cnt2; its_countmember(unique(yj(tmp)),yj(tmp))];
    Q(tmp,:)=[]; yj(tmp,:)=[];
end

%%% Entropy of A
p=cnt./n;
e=0;
for i=1:length(cnt)
   if p(i)== 0
   else
      e=e-p(i)*log(p(i));
   end   
end

%%% Entropy of B
p2=cnt2./n;
e2=0;
for i=1:length(cnt2)
   if p2(i)== 0
   else
      e2=e2-p2(i)*log(p2(i));
   end   
end

%%% Conditional Entropy
ce=e2-e;
end

function C = its_countmember(A,B)
% COUNTMEMBER - count members
%
%   C = COUNTMEMBER(A,B) counts the number of times the elements of array A are
%   present in array B, so that C(k) equals the number of occurences of
%   A(k) in B. A may contain non-unique elements. C will have the same size as A.
%   A and B should be of the same type, and can be cell array of strings.
%
%   Examples:
%     countmember([1 2 1 3],[1 2 2 2 2]) 
%        -> 1     4     1     0
%     countmember({'a','b','c'},{'a','x','a'}) 
%        -> 2     0     0
%
%   See also ISMEMBER, UNIQUE, HISTC

% for Matlab R13 and up
% version 1.2 (dec 2008)
% (c) Jos van der Geest
% email: jos@jasen.nl

% History:
% 1.0 (2005) created
% 1.1 (??): removed dum variable from [AU,dum,j] = unique(A(:)) to reduce
%    overhead
% 1.2 (dec 2008) - added comments, fixed some spelling and grammar
%    mistakes, after being selected as Pick of the Week (dec 2008)

error(nargchk(2,2,nargin)) ;

if ~isequal(class(A),class(B)),
    error('Both inputs should be the same class.') ;
end
if isempty(B),
    C = zeros(size(A)) ;
    return
elseif isempty(A),
    C = [] ;
    return
end

% which elements are unique in A, 
% also store the position to re-order later on
[AU,j,j] = unique(A(:)) ; 
% assign each element in B a number corresponding to the element of A
[L, L] = ismember(B,AU) ; 
% count these numbers
N = histc(L(:),1:length(AU)) ;
% re-order according to A, and reshape
C = reshape(N(j),size(A)) ;    
end
