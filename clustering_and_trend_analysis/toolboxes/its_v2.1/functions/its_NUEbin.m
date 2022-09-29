%% NON-UNIFORM EMBEDDING with BINNING

%%% INPUTS:
% data: N*M matrix of the M signals each having length N
% j: index (column) of the series considered as output, the one we want to describe
% candidates: two column vector of candidate series (col 1) and lag (col 2) indexes
% c: n. of quantization levels
% quantizza: set at 'n' to skip quantization (if you pass data already quantized)

% note: minshift>0 (typically 20) entails significance test using shift surrogates
% set minshift=0 to do random shuffling instead of time shift!

function out=its_NUEbin(data,j,candidates,c,numshift,minshift,alphashift,quantizza)

%% time series uniform quantization (UQ)
M=size(data,2);
if ~exist('quantizza'), quantizza='y'; end
if quantizza=='n'
    Yq=data;
else
    for m=1:M
        Yq(:,m)=its_quantization(data(:,m),c)-1;
    end
end

%% prepare for iteration
candidates=[candidates NaN*ones(size(candidates,1),1)]; %NaN in 3rd column

%Shannon Entropy of yj(n) based on UQ
Eyj=its_Ebin(Yq(:,j),c,'n');

%% compute Conditional entropy - iterate
exitcrit=0; %exit criterion: stay in loop if exitcrit==0
V=[];
ceOK(1)=Eyj; % la prima ce è l'entropia della target series
k=2;

while exitcrit==0
    disp(['NUembed, step ' int2str(k-1) '...']);
    ce=NaN*ones(size(candidates,1),1);
    
    % Test the candidates
    for i=1:size(candidates,1)
        if isnan(candidates(i,3))%test i-th candidate, only if not already included in V
            Vtmp=[V; candidates(i,1:2)];
            B=its_buildvectors(Yq,j,Vtmp);
            ce(i)=its_CEbin(B); %conditional entropy
        end
    end    
    
    %% select the candidate and update the embedding vector    
    ind_sel=find(ce==min(ce)); 
    if ~isempty(ind_sel)
        ind_sel=ind_sel(1); % index of the selected candidate
        candidates(ind_sel,3)=1; %mark as selected
        V=[V; candidates(ind_sel,1:2)]; %update embedding vector
        ceOK(k)=min(ce); %update the minimum ce at step k
    else    % ho testato tutti i candidati ed è sempre sceso!
        ceOK(k)=ceOK(k-1); %la fisso pianeggiante, per poter uscire
        V=[V; V(k-2,1:2)]; %fake update embedding vector (copio il precedente)
    end
    CMI(k-1)=ceOK(k-1)-ceOK(k); % conditional mutual information of the selected term
    
    %% test for significance of the new selected candidate (exit criterion)
    B=its_buildvectors(Yq,j,V); %prendo la observation matrix che ha dato la ce minima
    maxshift=size(B,1)-minshift;
    ces=NaN*ones(numshift,1);
    for is=1:numshift
        if minshift==0 % 26/3 alternativa a timshift: RANDOM SHUFFLING (minshift=0 as flag!)
            xs=B(:,size(B,2));
            xs=xs(randperm(length(xs)));
        else
            lagshift=fix(rand*(maxshift-minshift+1)+minshift);% shift casuale tra minshift e maxshift 
            xs=circshift(B(:,size(B,2)),lagshift);%is-esimo shift dei valori del candidato scelto       
        end             
        Bs=[B(:,1:size(B,2)-1) xs]; % sostituisco l'ultima colonna di B col termine shiftato
        ces(is)=its_CEbin(Bs);% conditional entropy per il surrogato       
    end
    CMIs=ceOK(k-1)-ces; %conditional mutual information of random new lag (according to surrogates)
    soglia(k-1)=prctile(CMIs,100*(1-alphashift));
    if CMI(k-1)<=soglia(k-1) % the term is not significant
        exitcrit=1;
    else % the term is significant: update k and stay in
        k=k+1;
        if ceOK(k-1)==0 % cannot be more predictable: exit
            ceOK(k)=0;
            exitcrit=1;
        end
    end
    
end %while

%%% After exiting: final values
ceOK=ceOK';
L=k-2; %optimal embedding length
VL=V(1:L,:); %optimal embedding

% CEmin=ceOK(k-1);

%%% arrange output structure
out.VL=VL;
out.Hy=Eyj;
out.Hy_VL=ceOK(L+1);
out.CE=ceOK;
out.CMI=CMI';
out.threshold=soglia';
out.VL1=V;


