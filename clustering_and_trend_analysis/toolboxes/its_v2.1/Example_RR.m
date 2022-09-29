%% Example of computation of univariate entropy measures on physiological time series
% Data: RR interval time series from a healthy subject (supine and upright)
% Measures: (Shannon) Entropy, Conditional Entropy, Information Storage
% Computation based on standard Uniform Embedding
% Estimators: model-based linear parametric, model-free based on binning, model-free based on kernels, model-free based on nearest neighbors
clear;close all;clc;

normalize='n'; % 'y' to study normalized time series (zero mean unit variance)

% parameters
p=3; % LINEAR: model order = number of past lagged components used for conditional entrpy analysis
c=6; % BINNING: n. of quantization levels
rfact=0.2; % KERNEL: threshold distance (fraction of SD of time series)
k=10; % NEAREST NEIGHBORS: number of neighbors
L=2; % BINNING, NEAREST NEIGHBORS: number of past lagged components used for conditional entrpy analysis (uniform embedding)

%% open heart period time series
data=load('upright.prn');
rr=data(:,1);

if normalize=='y'
    Y=(rr-mean(rr))./std(rr);
else
    Y=(rr-mean(rr)); % the mean should always be removed (this actually has effect only on the conditional entropy computed with linear estimator)
end


%% Estimation - linear method
% System Information
Hy1=its_Elin(Y); % Entropy

% New information
V1=its_SetLag(p,1);
B1=its_buildvectors(Y,1,V1); %observation matrix (for zero-mean data!)
Hy_y1=its_CElin(B1); % Conditional entropy

% Information storage
out1=its_SElin(Y,1,V1);
Sy1=out1.Sy;

%% Estimation - binning method
% System Information
Hy2=its_Ebin(Y,c); % Entropy

% New information
Yq=its_quantization(Y,c);
V2=its_SetLag(L,1); %indices of embedding vector
B2=its_buildvectors(Yq,1,V2); %observation matrix
Hy_y2=its_CEbin(B2); % Conditional entropy

% Information storage
out2=its_SEbin(Y,V2,1,c);
Sy2=out2.Sy;

%% Estimation - kernel method
% System Information
r=rfact*std(Y);
Hy3=its_Eker(Y,r); % Entropy

% New information
V3=its_SetLag(L,1); %indices of embedding vector
B3=its_buildvectors(Y,1,V3); %observation matrix
Hy_y3=its_CEker(B3,r); % Conditional entropy

% Information storage
out3=its_SEker(Y,V3,1,r);
Sy3=out3.Sy;

%% Estimation - nearest neighbor method
% System Information
Hy4=its_Eknn(Y,k); % Entropy

% New information
V4=its_SetLag(L,1); %indices of embedding vector
B4=its_buildvectors(Y,1,V4); %observation matrix
Hy_y4=its_CEknn(B4,k);

% Information storage
out4=its_SEknn(Y,V4,1,k);
Sy4=out4.Sy;



%% plots and disps
figure(1); plot(rr); title('RR time series');
xlabel('[beats]'); ylabel('[ms]');

clc;
disp(['ENTROPY:']);
disp(['Linear: Hy=' num2str(Hy1)]); 
disp(['Binning: Hy=' num2str(Hy2)]); 
disp(['Kernel: Hy=' num2str(Hy3)]); 
disp(['Nearest neighbor: Hy=' num2str(Hy4)]); 

disp(' ');
disp(['CONDITIONAL ENTROPY:']);
disp(['Linear: Hy|y=' num2str(Hy_y1)]); 
disp(['Binning: Hy|y=' num2str(Hy_y2)]); 
disp(['Kernel: Hy|y=' num2str(Hy_y3)]); 
disp(['Nearest neighbor: Hy|y=' num2str(Hy_y4)]); 

disp(' ');
disp(['INFORMATION STORAGE:']);
disp(['Linear: Sy=' num2str(Sy1)]); 
disp(['Binning: Sy=' num2str(Sy2)]); 
disp(['Kernel: Sy=' num2str(Sy3)]); 
disp(['Nearest neighbor: Sy=' num2str(Sy4)]); 





