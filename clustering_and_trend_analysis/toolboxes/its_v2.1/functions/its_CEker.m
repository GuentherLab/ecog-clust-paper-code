%% Kernel Estimation of the Conditional Entropy
% Computes the conditional information of the first column of B conditioned to the remaining columns (A)
%%% NOTE: this implementation is equivalent to the SAMPLE ENTROPY [Richman and Moorman]
% this uses maximum distance (Chebyshev) norm (Euclidian distance should also be implemented)
% r is threshold for distances, pass the absolute value

function Hy_y=its_CEker(B,r,norma)

if ~exist('norma','var') 
    norma='c'; %default Chebyshev
end

% if r not provided use rfact and std of data, otherwise use the r provided and don't care about rfact
% if ~exist('r','var') 
%     r=rfact*std(B(:,1));
% end

dataMat=flipdim(B',1);
dim=size(B,2)-1;
N=size(B,1)+dim;


% from here on, coded by Kijoon Lee,  kjlee@ntu.edu.sg [downloaded from MatlabCentral]
correl = zeros(1,2);
for m = dim:dim+1
    count = zeros(1,N-dim);
    tempMat = dataMat(1:m,:);
    
    for i = 1:N-m
        % calculate Chebyshev distance, excluding self-matching case
        if norma=='c'
            %%%dist = max(abs(tempMat(:,i+1:N-dim) - repmat(tempMat(:,i),1,N-dim-i)));
            dist = max(abs(tempMat(:,i+1:N-dim) - repmat(tempMat(:,i),1,N-dim-i)),[],1);
        else % Euclidian distance (check)
            tmp=tempMat(:,i+1:N-dim);
            dist=nan*ones(1,size(tmp,2));
            for q=1:size(tmp,2)
                dist(q)=norm(tempMat(:,i)-tmp(:,q),2);
            end
        end
        
        % calculate Heaviside function of the distance
        % User can change it to any other function
        % for modified sample entropy (mSampEn) calculation
        D = (dist < r);
        
        count(i) = sum(D)/(N-dim); %computes probability
    end
    
    correl(m-dim+1) = sum(count)/(N-dim); %average of probability
end


Hy_y = log(correl(1)/correl(2));



end


%%
% function saen = SampEnLee( dim, r, data, tau )
%%% dw da matlab  Kijoon Lee
% SAMPEN Sample Entropy
%   calculates the sample entropy of a given time series data

%   SampEn is conceptually similar to approximate entropy (ApEn), but has
%   following differences:
%       1) SampEn does not count self-matching. The possible trouble of
%       having log(0) is avoided by taking logarithm at the latest step.
%       2) SampEn does not depend on the datasize as much as ApEn does. The
%       comparison is shown in the graph that is uploaded.

%   dim     : embedded dimension
%   r       : tolerance (typically 0.2 * std)
%   data    : time-series data
%   tau     : delay time for downsampling (user can omit this, in which case
%             the default value is 1)
%
%---------------------------------------------------------------------
% coded by Kijoon Lee,  kjlee@ntu.edu.sg
% Mar 21, 2012
%---------------------------------------------------------------------


% if nargin < 4, tau = 1; end
% if tau > 1, data = downsample(data, tau); end
% clear; close all; clc
% data=(1:10)';
% r=0.2;
% dim=2;


% N = length(data);
% correl = zeros(1,2);
% dataMat = zeros(dim+1,N-dim);
% for i = 1:dim+1
%     dataMat(i,:) = data(i:N-dim+i-1);
% end
% 
% for m = dim:dim+1
%     count = zeros(1,N-dim);
%     tempMat = dataMat(1:m,:);
%     
%     for i = 1:N-m
%         % calculate Chebyshev distance, excluding self-matching case
%         dist = max(abs(tempMat(:,i+1:N-dim) - repmat(tempMat(:,i),1,N-dim-i)));
%         
%         % calculate Heaviside function of the distance
%         % User can change it to any other function
%         % for modified sample entropy (mSampEn) calculation
%         D = (dist < r);
%         
%         count(i) = sum(D)/(N-dim);
%     end
%     
%     correl(m-dim+1) = sum(count)/(N-dim);
% end
% 
% 
% saen = log(correl(1)/correl(2));
% end