%% SPECTRAL AND TRANSFER FUNCTION MATRICES OF VAR PROCESSES

%%% inputs: 
% Am=[A(1)...A(p)]: M*pM matrix of the MVAR model coefficients (strictly causal model)
% Su: M*M covariance matrix of the input noises
% N= number of points for calculation of the spectral functions (nfft)
% Fs= sampling frequency

%%% outputs:
% S= Spectral Matrix
% H= Tranfer Function Matrix
% f= frequency vector

function [S,H,f] = var_spectra(Am,Su,N,Fs)

M= size(Am,1); % Am has dim M*pM
p = size(Am,2)/M; % p is the order of the VAR model

if nargin<2, Su = eye(M,M); end; % if not specified, we assume uncorrelated noises with unit variance as inputs 
if nargin<3, N = 512; end;
if nargin<4, Fs= 1; end;     
if all(size(N)==1),	 %if N is scalar
    f = (0:N-1)*(Fs/(2*N)); % frequency axis
else            % if N is a vector, we assume that it is the vector of the frequencies
    f = N; N = length(N);
end;

s = exp(i*2*pi*f/Fs); % vector of complex exponentials
z = i*2*pi/Fs;


%% Initializations: spectral matrices have M rows, M columns and are calculated at each of the N frequencies
S=zeros(M,M,N); % Spectral Matrix

A = [eye(M) -Am]; % matrix from which M*M blocks are selected to calculate spectral functions
invSu=inv(Su);

%% computation of spectral functions
for n=1:N, % at each frequency

        As = zeros(M,M); % matrix As(z)=I-sum(A(k))
        for k = 1:p+1,
            As = As + A(:,k*M+(1-M:0))*exp(-z*(k-1)*f(n));  %indicization (:,k*M+(1-M:0)) extracts the k-th M*M block from the matrix B (A(1) is in the second block, and so on)
        end;
        
        %%% Transfer matrix
        H(:,:,n)  = inv(As);
        
        %%% Spectral matrix
        S(:,:,n)  = H(:,:,n)*Su*H(:,:,n)'; % ' stands for Hermitian transpose
        
end;


