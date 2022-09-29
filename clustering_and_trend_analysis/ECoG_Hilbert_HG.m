function ECoG = ECoG_Hilbert_HG(ECoG, processing)

% band pass + calcualte hilbert and average
fprintf('---- High gamm calculation  \n')

switch processing
    case 'highgamma'
        edges = logspace(log10(70), log10(150), 9);
    case 'beta'
        edges = logspace(log10(15), log10(30), 6);
end

Ast1 = 60;
Ap = 1;
Ast2 = 60;
Fs = ECoG.srate;
for i = 1:length(edges)-1

    Fst1 = edges(i) - 1;
    Fp1 = edges(i);
    Fp2 = edges(i+1);
    Fst2 = edges(i+1) + 1;
    
    H = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2,Fs);
    Hd{i} = design(H,'butter');
    
%     [h{i}, w{i}] = freqz(Hd{i});
%     plot(w{i}/pi * Fs/2, 20*log10(abs(h{i})));
    
end
% Hd = design(H,'equiripple');
% fvtool(Hd);
for ch =1 :size(ECoG.data,1)
    
    tempSig = [zeros(1,10000) ECoG.data(ch,:) zeros(1,10000)];
    hilbert_tempVar = [];
    for i = 1:length(edges)-1
        
        tempVar= filtfilt(Hd{i}.sosMatrix,Hd{i}.ScaleValues,tempSig);
        tempVar1 = tempVar(10001:end-10000);
        hilbert_tempVar(i,:) = abs(hilbert(tempVar1));
    end
    ECoG.data_hilbert(ch,:) = mean(hilbert_tempVar,1);
end


% Ayoub's orignal way, band edges not correct
% 
% fprintf('---- High gamm calculation  \n')
% CentralFreq = logspace(log10(70),log10(150),8);
% bw = diff(CentralFreq)/2;
% lf = [bw(1) bw];
% hf = [bw bw(end)];
% 
% for i = 1:8
%     Fst1 = CentralFreq(i)-lf(i)-1;
%     Fp1 = CentralFreq(i)-lf(i);
%     Fp2 = CentralFreq(i)+hf(i);
%     Fst2 = CentralFreq(i)+hf(i)+1;
%     Ast1 = 60;
%     Ap = 1;
%     Ast2 = 60;
%     Fs = 1000;
%     H = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2,Fs);
%     Hd{i} = design(H,'butter');
% end
% % Hd = design(H,'equiripple');
% % fvtool(Hd);
% for ch =1 :size(ECoG.data,1)
%     
%     tempSig = [zeros(1,10000) ECoG.data(ch,:) zeros(1,10000)];
%     for i = 1:8
%         
%         tempVar= filtfilt(Hd{i}.sosMatrix,Hd{i}.ScaleValues,tempSig);
%         tempVar1 = tempVar(10001:end-10000);
%         hilbert_tempVar(i,:) = abs(hilbert(tempVar1));
%     end
%     ECoG.data_hilbert(ch,:) = mean(hilbert_tempVar,1);
% end
