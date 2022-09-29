clear;
close all;
clc;


sub_str = 'S376'; % Need Separate file for S376
% Changed from Ayoubs to get functionality here, mostly just going off
% trigger

new_Fs = 1000;

base_path = error('Enter data directory path here');

subject_path =  '/projectnb/busplab/UserData/adaliri/RepSup/ProcessedDataset/';
% subject_path =  '/projectnb/busplab/UserData/adaliri/RepSup/NewProcessedDataset';

load(fullfile(base_path, sub_str, 'ECoGALL.mat'));

stimuliOnset = struct([]);
for k = 1:length(ECoGALL)

    % Get local data copy
    ECoG  = ECoGALL{k};
    
    % Trigger sampling frequency
    Fs = ECoG.Fs(end)
    
    % Calculate Trigger
    sigData = double(ECoG.Trigger) - double(ECoG.Speaker); 
    sigData= sigData - mean(sigData);
    %         sigData(1:12*Fs) = 0;
    tempSig = sigData;
    winConv = ones(1, 5*round(Fs/new_Fs));
    envData = [];
    envData = conv(abs(tempSig),winConv);
%     envData = tempSig;

    thrVal = 2000;
    tDurMin = .05;
    envData(envData<thrVal)= 0;
    envData(envData>=thrVal) = 1;
    voiceSig = diff(envData);
    aInd = find(voiceSig>0);
    durMin = round(tDurMin*Fs);
    kk = 1;
    trigInd = [];
    for i = 1:length(aInd)
        if (aInd(i)+durMin) < length(envData)
            if sum(envData(aInd(i):aInd(i)+durMin )) == durMin
                trigInd(kk) = aInd(i);
                kk = kk +1;
            end
        end
    end



    figure('Name',sub_str,'Position',[10 10 1800 900])
    to = linspace(0,length(sigData)/ECoG.Fs(end),length(sigData));
    plot(to,sigData/max(abs(sigData)),'color',[1 1 1]/3);
    hold on;
    stem(to(trigInd),.4*ones(size(to(trigInd))),'ro','LineWidth',3)






    % calculate onset time
    
    sigData = double(ECoG.Mic);
    sigData(1:5*Fs) = mean(sigData);
    sigData(350*Fs:end) = mean(sigData);
    sigData= sigData - mean(sigData);

    %         [b1,a1] = butter(8,100/(Fs/2),'high');
    %         [b2,a2] = butter(8,2000/(Fs/2),'low');
    %
    %         tempSig = filtfilt(b1,a1,sigData);
    %         tempSig = filtfilt(b2,a2,tempSig);
    tempSig = sigData;
    winConv = ones(1,10*round(Fs/new_Fs));
    envData = [];
    envData = conv(abs(tempSig),winConv);



    thrVal = 2000;
    tDurMin = .3;
    envData(envData<thrVal)= 0;
    envData(envData>=thrVal) = 1;
    voiceSig = diff(envData);
    aInd = find(voiceSig>0);
    durMin = round(tDurMin*Fs);
    kk = 1;
    onsetInd = [];
    for i = 1:length(aInd)
        if (aInd(i)+durMin) < length(envData)
            if sum(envData(aInd(i):aInd(i)+durMin )) == durMin
                onsetInd(kk) = aInd(i);
                kk = kk +1;
            end
        end
    end



    figure('Name',sub_str,'Position',[10 10 1800 900])
    to = linspace(0,length(sigData)/ECoG.Fs(end),length(sigData));
    plot(to,sigData/max(abs(sigData)),'color',[1 1 1]/3);
    hold on;
    stem(to(onsetInd),.4*ones(size(to(onsetInd))),'ro','LineWidth',3)




    % calculate offset time
    tempSig = tempSig(end:-1:1);
    envData = conv(winConv,abs(tempSig));

    thrVal = .005;
    tDurMin = .2;

    envData(envData<thrVal)= 0;
    envData(envData>=thrVal) = 1;
    voiceSig = diff(envData);
    aInd = find(voiceSig>0);
    durMin = round(tDurMin*Fs);
    kk = 1;
    offInd = [];
    for i = 1:length(aInd)
        if (aInd(i)+durMin) < length(envData)
            if sum(envData(aInd(i):aInd(i)+durMin )) == durMin
                offInd(kk) = aInd(i);
                kk = kk +1;
            end
        end
    end

    offInd = length(tempSig) - offInd;
    stem(to(offInd),.4*ones(size(to(offInd))),'go','LineWidth',3)
    stem(to(trigInd),.4*ones(size(to(trigInd))),'ko','LineWidth',3)

    % onset times for jeremy
    stimuliOnset(k).VisualOnset = trigInd;
    stimuliOnset(k).SpeechOnset = onsetInd;

%     % converting to ms (new_Fs = 1000)
%     trigInd  = round(new_Fs * trigInd/Fs);
%     onsetInd = round(new_Fs * onsetInd/Fs);
%     offInd   = sort(round(new_Fs * offInd/Fs));
% 
%     % showing on the command window for inspection purposes
%     [numel(trigInd) numel(onsetInd)  numel(offInd) ]
% 
%     % save the values 
%     ECoG.trigInd = trigInd;
%     ECoG.onsetTime = onsetInd;
%     %         EEG.offsetTime = offInd;
% 
%     % careate a table for pairs
%     A = trigInd(:);
%     onsetTable(:,1:2) = reshape(A,2,length(A)/2)';
%     A = onsetInd(:);
%     onsetTable(:,3:4) = reshape(A,2,length(A)/2)';
%     %         A = offInd(:);
%     %         onsetTable(:,5:6) = reshape(A,2,36)';
%     ECoG.onsetTable = onsetTable;
% 
%     ECoGALL{k} = ECoG;
%     ECoG = [];


end

% %     outName = fullfile(subject_path,BlockFileName);
% % Save off updated ECoG
% outName = fullfile(base_path, sub_str, 'ECoGALL.mat');
% save(outName,'ECoGALL','-v7.3');

% Save off onset info
outName = fullfile(base_path, sub_str, [sub_str 'OnsetTime.mat']);
save(outName,'stimuliOnset');






