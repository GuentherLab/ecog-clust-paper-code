%% THEORETICAL COEFFICIENTS FOR SIMULATED MVAR PROCESSES

function [Am,Su,Ak]=var_simulations(M,par)
% clear; close all;clc;
% M=2;
% par.poles=([0.8 0.1; 0.9 0.3]); % Autonomous oscillations in each process
% par.coup=[1 2 1 0.5; 2 1 2 -0.2]; % in each row: "i j k c" to impose coupling from i to j at lag k with coeff c 
% par.Su=[1 1]; % variance of innovation processes

if isempty(par.coup)
    p=2;
else
    p=max(2,max(par.coup(:,3)));
end

Su=eye(M); %innovation covariance matrix (diagonal)
for k=1:M
    Su(k,k)=par.Su(k);
end

Ak=zeros(M,M,p);
for k=1: size(par.poles,1)
    Ak(k,k,1)=2*par.poles(k,1)*cos(2*pi*par.poles(k,2));
    Ak(k,k,2)=-par.poles(k,1)^2;
end
for k=1:size(par.coup,1)
    Ak(par.coup(k,2),par.coup(k,1),par.coup(k,3))=par.coup(k,4);
end

Am=[];
for kk=1:p
    Am=[Am Ak(:,:,kk)];
end

% stability check
E=eye(M*p);AA=[Am;E(1:end-M,:)];lambda=eig(AA);lambdamax=max(abs(lambda));
if lambdamax>=1
    error('The simulated VAR process is not stable');
end

end


