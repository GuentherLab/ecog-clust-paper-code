%% SELF ENTROPY (INFORMATION STORAGE) for LINEAR GAUSSIAN PROCESSES (UNIFORM EMBEDDING)
% computes the self-entropy of a time series (jj-th column of data matrix)

%%% INPUTS:
% data: N*M matrix of M time series each having length N
% jj: index of target variable - the one of which we want to compute the SE
% V: vector of candidates (to be pre-determined from uniform embedding)

function out=its_SElin(data,jj,V)

for m=1:size(data,2)
    data(:,m) = data(:,m) - mean(data(:,m)); % entropy linear is sensitive to mean - we remove it
end

% Entropy
[Hy, varY]=its_Elin(data(:,jj));
% Hy=ret.Hy;
% varY=ret.SigmaY;

if isempty(V) % no conditioning: SE=0, CE=E
    [out.Hy_y,out.Sy_y]=its_Elin(data(:,jj));
    out.Hy=Hy;
    out.Hy_y=Hy;
    out.SigmaY=varY;
    out.SigmaY_y=varY;
    out.Sy=0;
    return
end



% form the observation matrix
B=its_buildvectors(data,jj,V);
A=B(:,2:end);
tmp=V(:,1); 
i_Y= tmp==jj;
M_y=B(:,1);
M_Y=A(:,i_Y);
M_yY=[M_y M_Y];

%%% Regression
[Hy_y,Sy_y,Uy_y,Am]=its_CElin(M_yY);

% Self Entropy
SE=Hy-Hy_y;


Uy=data(:,jj); % residuals of unconditioned regression are the series itself
pval = its_LinReg_Ftest(Uy_y,Uy,size(M_Y,2),0);

%%% output structure
out.Sy=SE;
out.p_Sy=pval;
out.Hy=Hy;
out.Hy_y=Hy_y;
out.Uy_y=Uy_y;
out.coeff=Am;
out.SigmaY=varY;
out.SigmaY_y=Sy_y;

end




