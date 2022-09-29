function ECoG = Early_Preprocessing(ECoG, line_noise_exit)

if nargin < 2
    line_noise_exit = false;
end

addpath(genpath('EEG-Clean-Tools-master'));
addpath(genpath('eeglab13_5_4b'));

%% load data
% pathDATA = '/projectnb/busplab/UserData/adaliri/ecogData/dataset/';
% load (fullfile(ECoG_File));
% fprintf('Starting to process subject : %s  \n',ECoG_File)
% High pass filter + Detrend the signal for detecting bad channels


fprintf('---- High pass filter + Detrend \n')
% Fst1 = .5;
% Fp1 = 3;
% Fp2 = 170;
% Fst2 = 175;
% Ast1 = 60;
% Ap = 1;
% Ast2 = 60;
% Fs = 1000;
% H = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2,Fs);
% Hd = design(H,'butter');
% % fvtool(Hd);
% for ch = 1:size(ECoG.data,1)
%     tempSig = [zeros(1,100000) ECoG.data(ch,:) zeros(1,100000)];
%     tempVar= filtfilt(Hd.sosMatrix,Hd.ScaleValues,tempSig);
%     ECoG.data(ch,:) = tempVar(100001:end-100000);
% end
ECoG.data = ECoG.data - repmat(mean(ECoG.data,2),1,size(ECoG.data,2));


%% LineNoise Removal (PREP)
fprintf('---- Line Noise removal \n')
ECoG.srate = 1000;
ECoG.trials = 1;
ECoG.nchan=size(ECoG.data,1);
[ECoG, lineNoiseOut] = cleanLineNoise(ECoG);

if line_noise_exit
    ECoG.badChan = [];
    return
end
    
%% Bad Cahennels removal
fprintf('---- Bad channel detection  \n')
kurtThre = 8;
kurtVal = kurtosis(ECoG.data');

goodChan = find(kurtVal<kurtThre);
ECoG.badChan = find(kurtVal>=kurtThre);

% fprintf('---- Bad channel detection  \n')
% kurtThre = 3;
% tempSig = (EEG.data-repmat(trimmean(EEG.data,20,2),1,size(EEG.data,2)))./repmat(std(EEG.data,[],2),1,size(EEG.data,2));
% [EEG.Kurt, EEG.badChan] = rejkurt( tempSig, kurtThre, [], 2);


%% estimate TRUE reference and re-reference
% refCAR = mean(ECoG.data(goodChan,:));
fprintf('---- Re-referencing  \n')
refCAR = trimmean(ECoG.data,20);
ECoG.data = ECoG.data - repmat(refCAR,size(ECoG.data,1),1);





end


