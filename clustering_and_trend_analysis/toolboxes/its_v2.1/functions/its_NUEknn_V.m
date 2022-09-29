%% NONUNIFORM MIXED EMBEDDING for PREDICTIVE INFORMATION based on nearest neighbors
%% Approach 2-mod): Conditional Mutual Information estimation
%% VERSION WITH MULTIVARIATE (VECTOR) TARGET: j can be vector variable
% w.r.t. MInu_knn, this is modified to select components through CMI, instead MI
% stops at minimum CMI (like in Kugiumtzis2013)

%%% INPUTS:
% data: N*M matrix of the M signals each having length N
% j: index (column) of the series considered as output, the one we want to describe
% candidates: two column vector of candidate series (col 1) and lag (col 2) indexes
% k: number of the nearest neighbor to search for distance
% metric: 'maximum' or 'euclidian'

% note: minshift>0 (typically 20) entails significance test using shift surrogates
% set minshift=0 to do random shuffling instead of time shift!

function out=its_NUEknn_V(data,j,candidates,nnk,num_rnd,minshift,alpha_rnd,metric)

if ~exist('metric','var'), metric='maximum'; end

%% prepare for iteration
candidates=[candidates NaN*ones(size(candidates,1),1)]; %NaN in 3rd column
dimj=length(j);
%Shannon Entropy of yj(n) based on knn
% Eyj=SEknn(data(:,j),nnk,metric);


%% compute CMI - iterate
exitcrit=0; %exit criterion: stay in loop if exitcrit==0
V=[];
cmiOK(1)=0; % initialize MI between target series and embedding vector
k=2;

while exitcrit==0
    disp(['NUembed, step ' int2str(k-1) '...']);
    c_mi=NaN*ones(size(candidates,1),1);
    %d=NaN*ones(size(candidates,1),1);
    
    % Test the candidates
    for i=1:size(candidates,1)
        if isnan(candidates(i,3)) % test i-th candidate, only if not already included in V
            Vtmp=[V; candidates(i,1:2)];
            B=its_buildvectors(data,j,Vtmp);
            c_mi(i)=its_CMIknn_V(B,dimj,nnk,metric); % CONDITIONAL mutual information knn
        end
    end    
    
    %% select the candidate and update the embedding vector    
    ind_sel=find(c_mi==max(c_mi)); 
    if ~isempty(ind_sel)
        ind_sel=ind_sel(1); % index of the selected candidate
        candidates(ind_sel,3)=1; %mark as selected
        V=[V; candidates(ind_sel,1:2)]; %update embedding vector
        cmiOK(k)=max(c_mi); %update the maximum c_mi at step k
    else    % ho testato tutti i candidati ed è sempre salito!
        cmiOK(k)=cmiOK(k-1); %la fisso pianeggiante, per poter uscire
        V=[V; V(k-2,1:2)]; %fake update embedding vector (copio il precedente)
    end
    
    %% test for significance of the new selected candidate (exit criterion)
    B=its_buildvectors(data,j,V); %riprendo la observation matrix che ha dato la c_mi massima
    CMI(k-1)=its_CMIknn_V(B,dimj,nnk,metric); % conditional mutual information of the selected term (nn estimate)

    maxshift=size(B,1)-minshift;
    CMIs=NaN*ones(num_rnd,1);
    for is=1:num_rnd
        if minshift==0 % alternative to timeshift: RANDOM SHUFFLING (minshift=0 used as flag)
            xs=B(:,size(B,2));
            xs=xs(randperm(length(xs)));
            clear x1s;
            for jjs=1:dimj
                x1s1=B(:,jjs);
                x1s(:,jjs)=x1s1(randperm(length(x1s1))); % Dimitris suggestion: randomize also target variable
            end
        else
            lagshift=fix(rand*(maxshift-minshift+1)+minshift);% shift casuale tra minshift e maxshift 
            xs=circshift(B(:,size(B,2)),lagshift);%is-esimo shift dei valori del candidato scelto
            lagshift1=fix(rand*(maxshift-minshift+1)+minshift);
            clear x1s;
            for jjs=1:dimj
                x1s(:,jjs)=circshift(B(:,jjs),lagshift1); % Dimitris suggestion: randomize also target variable
            end
        end             
        Bs=B; Bs(:,end)=xs; % sostituisco l'ultima colonna di B col termine shiftato
        Bs(:,1:dimj)=x1s; % Dimitris suggestion: randomize also target variable
        CMIs(is)=its_CMIknn_V(Bs,dimj,nnk,metric);
    end
    soglia(k-1)=prctile(CMIs,100*(1-alpha_rnd));
    
    if CMI(k-1)<=soglia(k-1) % the term is not significant
        exitcrit=1;
    else % the term is significant: update k and stay in
        k=k+1;
    end
        
    
end %while

%%% After exiting: final values
cmiOK=cmiOK';
L=k-2; %optimal embedding length
VL=V(1:L,:); %optimal embedding

%%% arrange output structure
out.VL=VL;
out.CMIok=cmiOK';
out.CMI=CMI';
out.threshold=soglia';
out.VL1=V;


end



