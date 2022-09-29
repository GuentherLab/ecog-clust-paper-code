%% Compares linear, knn, binning and kernels for estimating entropy, conditional entropy and information storage of a Gaussian AR process
clear; clc;

%% Parameters
numfig=1; % figure where to plot
N=1000;% length of simulated time series
which_estimator='bin'; % 'lin', 'bin', 'ker', 'knn'

kv=[5 10 30]; % for knn estimator, n. of neighbors
qv=[3 6 12]; % for bin estimator,n. of bins
rfactv=[0.1 0.2 0.5]; % for ker estimator, percent of SD for kernel method

numreal=10;% number of ralizations of the simulations to be generated

%%% Simulated AR process
rhov=[0 0.4 0.6 0.8 0.9]; % pole modulus
f1=0.25; % AR process for Y1

V=[1 1; 1 2]; %embedding vector (the optimal for this simu: 2 points in the past)


%% computations
numrho=length(rhov);

if which_estimator(1)=='l', numpar=1; else numpar=length(kv); end
eE=nan*ones(numrho,numreal,numpar); eCE=eE; eSE=eE;
for is=1:numrho
    
    rho=rhov(is);
    %%% simulated AR(2) process (one pole)
    Am=[2*rho*cos(2*pi*f1) -rho^2];
    Su=1;

    % theoretical information dynamics
    numlags=2; reto = its_CElinVAR1(Am,Su,numlags);
    H(is)=0.5*log(2*pi*exp(1)); % since I will work on normalized data, entropy is the same for any variance
    CE(is)=0.5*log(2*pi*exp(1)*reto.Sy_y/reto.Sy); 
    SE(is)=H(is)-CE(is);  
    
    % Realizations
    for ir=1:numreal
        clc; disp(['configuration ' int2str(is) ', realization ' int2str(ir)]);
        
        %%% generate data
        E=sqrt(Su)*randn(N,1);
        Yo=var_filter(Am,E);
        
        % make the simulated signal of zero mean and unit variance
        x = (Yo-mean(Yo)) ./ std(Yo);
        
        
        %%% ANALYSIS
        B=its_buildvectors(x,1,V); %Observation Matrix       
        switch which_estimator
            case 'lin'
                eE(is,ir)=its_Elin(x);
                eCE(is,ir)=its_CElin(B);
                out=its_SElin(x,1,V);
                eSE(is,ir)=out.Sy;
                
            case 'bin'
                for ip=1:numpar
                    q=qv(ip);
                    eE(is,ir,ip)=its_Ebin(x,q);
                    xq=its_quantization(x,q)-1; Bq=its_buildvectors(xq,1,V); %Observation Matrix
                    eCE(is,ir,ip)=its_CEbin(Bq);
                    out=its_SEbin(x,V,1,q);
                    eSE(is,ir,ip)=out.Sy;
                end
                
            case 'ker'
                for ip=1:numpar
                    rfact=rfactv(ip);
                    r=rfact*std(x);
                    eE(is,ir,ip)=its_Eker(x,r);
                    eCE(is,ir,ip)=its_CEker(B,r);
                    out=its_SEker(x,V,1,r);
                    eSE(is,ir,ip)=out.Sy;
                    
                end
                
            case 'knn'
                 for ip=1:numpar
                    k=kv(ip);
                    eE(is,ir,ip)=its_Eknn(x,k,'maximum');
                    eCE(is,ir,ip)=its_CEknn(B,k,'maximum');
                    outk=its_SEknn(x,V,1,k,'maximum');
                    eSE(is,ir,ip)=outk.Sy;
                end
        end
        
    end

end

%% theoretical (exact) values of information dynamics
rhovt=(0:0.01:0.9);
for it=1:length(rhovt)
    rho=rhovt(it);
    %%% simulated AR(2) process (one pole)
    Am=[2*rho*cos(2*pi*f1) -rho^2];
    Su=1;
    % theoretical information dynamics
    numlags=2; reto = its_CElinVAR1(Am,Su,numlags);
    Ht(it)=0.5*log(2*pi*exp(1)); % since I will work on normalized data, entropy is the same for any variance
    CEt(it)=0.5*log(2*pi*exp(1)*reto.Sy_y/reto.Sy); 
    SEt(it)=Ht(it)-CEt(it);
end


%% figure
%%% distributions
lo=10; hi=90;

eCE_prc=prctile(eCE,[lo 50 hi],2);
eSE_prc=prctile(eSE,[lo 50 hi],2);

 label_legend{1}='true';
switch which_estimator
    case 'lin'
        eE_prc=prctile(eE,[lo 50 hi],2);
        label_legend{2}='estimated';

    case 'bin'
        for ip=1:numpar
            label_legend{ip+1}=['q=' int2str(qv(ip))];
            eE_prc(:,:,ip)=prctile(eE(:,:,ip),[lo 50 hi],2);
            eCE_prc(:,:,ip)=prctile(eCE(:,:,ip),[lo 50 hi],2);
            eSE_prc(:,:,ip)=prctile(eSE(:,:,ip),[lo 50 hi],2);
        end

    case 'ker'
        for ip=1:numpar
            label_legend{ip+1}=['r%=' num2str(rfactv(ip))];
            eE_prc(:,:,ip)=prctile(eE(:,:,ip),[lo 50 hi],2);
            eCE_prc(:,:,ip)=prctile(eCE(:,:,ip),[lo 50 hi],2);
            eSE_prc(:,:,ip)=prctile(eSE(:,:,ip),[lo 50 hi],2);
        end


    case 'knn'
        for ip=1:numpar
            label_legend{ip+1}=['k=' int2str(kv(ip))];
            eE_prc(:,:,ip)=prctile(eE(:,:,ip),[lo 50 hi],2);
            eCE_prc(:,:,ip)=prctile(eCE(:,:,ip),[lo 50 hi],2);
            eSE_prc(:,:,ip)=prctile(eSE(:,:,ip),[lo 50 hi],2);
        end
            
end


colori='krbcm';
dx=0.02;
figure(numfig); clf;
set(gcf, 'Position', [1 1 1000 800]);
subplot(1,3,1);
plot(rhovt,Ht,'g-'); hold on;
for ip=1:numpar
    errorbar(rhov+(ip-ceil(numpar/2))*dx,eE_prc(:,2,ip),eE_prc(:,2,ip)-eE_prc(:,1,ip),eE_prc(:,3,ip)-eE_prc(:,2,ip), [colori(ip) '.']);
end
xlim([min(rhovt)-2*dx max(rhovt)+2*dx]);
set(gca,'XTick',rhov); set(gca,'XTickLabel',rhov);
title(['H, ' which_estimator]);
legend(label_legend,'location','best');
ylim([min([min(min((eE_prc(:,1,:)))) 0]) max([max(max((eE_prc(:,3,:)))) 1.5])]);
xlabel('\rho')

subplot(1,3,2);
plot(rhovt,CEt,'g-'); hold on;
for ip=1:numpar
errorbar(rhov+(ip-ceil(numpar/2))*dx,eCE_prc(:,2,ip),eCE_prc(:,2,ip)-eCE_prc(:,1,ip),eCE_prc(:,3,ip)-eCE_prc(:,2,ip), [colori(ip) '.']);
end
xlim([min(rhovt)-2*dx max(rhovt)+2*dx]);
set(gca,'XTick',rhov); set(gca,'XTickLabel',rhov);
title(['CE, ' which_estimator]);
legend(label_legend,'location','best');
ylim([min([min(min((eCE_prc(:,1,:)))) 0]) max([max(max((eCE_prc(:,3,:)))) 1.5])]);
xlabel('\rho')

subplot(1,3,3);
plot(rhovt,SEt,'g-'); hold on;
for ip=1:numpar
errorbar(rhov+(ip-ceil(numpar/2))*dx,eSE_prc(:,2,ip),eSE_prc(:,2,ip)-eSE_prc(:,1,ip),eSE_prc(:,3,ip)-eSE_prc(:,2,ip), [colori(ip) '.']);
end
xlim([min(rhovt)-2*dx max(rhovt)+2*dx]);
set(gca,'XTick',rhov); set(gca,'XTickLabel',rhov);
title(['IS, ' which_estimator]);
legend(label_legend,'location','best');
ylim([min([min(min((eSE_prc(:,1,:)))) 0]) max([max(max((eSE_prc(:,3,:)))) 0.75])]);
xlabel('\rho')
