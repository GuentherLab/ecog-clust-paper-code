%% Compare uniform and non-uniform embedding on simulated coupled henon maps
clear;close all;clc;

%% parameters
which_estimator='knn'; % 'ker', 'knn'
k=10; % if, knn, n. of neighbors
rfact=0.2; % if kernel, threshold r
Lmax=5; %number of lags for nonuniform embedding scheme

num_rnd=100; minshift=0;alpha_rnd=0.05; % CMI-based NUembed: surro-based termination criterion 

% single_analysis=0; %flag to execute only one PTE, or all together
% full_analysis=1-single_analysis;

%% Simulation - Multiple Henon [Kugiumtzis, PRE 2013]
M=5; % number of systems
C=0.5; %coupling strength
N=500; %series length
Yo = sim_coupledhenonmaps2(M,C,N);


%% Non Uniform Embedding
%%% Normalization - to be always done for knn
for m=1:M 
    Y(:,m)=(Yo(:,m)-mean(Yo(:,m)))/std(Yo(:,m)); 
end

% set candidates
pV=Lmax*ones(1,M);
tau=ones(1,M);
candidates=its_SetLag(pV,tau);
VLue=candidates; % uniform embedding vector


%% full analysis
PTEue=NaN*ones(M,M); % Partial Transfer Entropy, uniform embedding
PTEnue=NaN*ones(M,M); % Partial Transfer Entropy, nonuniform embedding
for jj=1:M
    clc; disp([which_estimator ', target ' int2str(jj) ' of ' int2str(M) ', embedding...']);
    
    %nonuniform embedding of target time series
    switch which_estimator
        case 'ker'
            out=its_NUEker(Y,jj,candidates,rfact*std(Y(:,jj)),num_rnd,minshift,alpha_rnd);
        case 'knn'
            out=its_NUEknn(Y,jj,candidates,k,num_rnd,minshift,alpha_rnd); 
    end  
    VLnue{jj}=out.VL; % nonuniform embedding vector
    
    for ii=1:M
        if ii~=jj
            switch which_estimator
                case 'ker'
                    outTEue=its_PTEker(Y,VLue,ii,jj,rfact*std(Y(:,jj)));
                    outTEnue=its_PTEker(Y,VLnue{jj},ii,jj,rfact*std(Y(:,jj)));
                case 'knn'
                    outTEue=its_PTEknn(Y,VLue,ii,jj,k);
                    outTEnue=its_PTEknn(Y,VLnue{jj},ii,jj,k); 
            end  
            PTEue(jj,ii)=outTEue.Txy_z;
            PTEnue(jj,ii)=outTEnue.Txy_z;
        end
    end

end

%% Plot results
figure(1); clf;
for m=1:M
    subplot(M,1,m); plot(Y(:,m)); zoom xon;
end

ymin=min([min(PTEue(:)) min(PTEnue(:))]);
ymax=max([max(PTEue(:)) max(PTEnue(:))]);

figure(2); clf;
subplot(1,2,1);
imagesc(PTEue, [ymin ymax]); xlabel('i'); ylabel('j'); title('TE_{i \rightarrow j | (1,...,M)\\i}, uniform embedding');
colorbar
set(gca,'Xtick',1:5); set(gca,'Ytick',1:5);

subplot(1,2,2);
imagesc(PTEnue, [ymin ymax]); xlabel('i'); ylabel('j'); title('TE_{i \rightarrow j | (1,...,M)\\i}, NONuniform embedding');
colorbar
set(gca,'Xtick',1:5); set(gca,'Ytick',1:5);

set(gcf, 'Position', [1 1 1000 400]);

