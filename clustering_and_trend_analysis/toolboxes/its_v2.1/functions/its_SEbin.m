%% SELF ENTROPY WITH BINNING and NON-UNIFORM EMBEDDING
% Computes the self entropy of a scalar target variable (column jj of Y)
% V= assigned embedding vector (only the components of jj are used)

%%% inputs are:
% Y: data matrix (N*M)
% jj: index of target variable 
% V: vector of candidates (to be pre-determined from uniform embedding or using CEnu_bin.m)
% c: n. of quantization levels for binning
% quantizza: set at 'n' to skip quantization (if you pass data already quantized)

function out=its_SEbin(Y,V,jj,c,quantizza)

Y=Y(:,jj); %we will use only the target

% time series uniform quantization (UQ)
if ~exist('quantizza'), quantizza='y'; end
if quantizza=='n'
    Yq=Y;
else
    Yq=its_quantization(Y,c)-1;
end

Hy=its_Ebin(Yq(:,jj),c,'n');

if isempty(V)
    out.Hy=Hy;
    out.Hy_y=Hy;
    out.Sy=0;
%     SE=0;
    return
end


tmp=V(:,1);
V=V(tmp==jj,:); %remove components other than jj

M_yY=its_buildvectors(Yq,1,V);
M_Y=M_yY(:,2:end);

if isempty(M_Y)
    Hy_y=Hy;
else
    Hy_y=its_CEbin(M_yY);
end

Sy=Hy-Hy_y;

out.Hy=Hy;
out.Hy_y=Hy_y;
out.Sy=Sy;


end










