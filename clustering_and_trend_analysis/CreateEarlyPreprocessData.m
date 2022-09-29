clear all;
close all;
clc;

addpath('PreprocessTools');
addpath(genpath('EEG-Clean-Tools-master'));
addpath(genpath('eeglab13_5_4b'));
    
sub_str = 'S376'; %{'S352', 'S357', 'S362', 'S369', 'S372', 'S376'}
line_noise_exit_flag = true;

basePath = error('Enter data directory path here');
dataPath = fullfile(basePath, 'RawData');
subPath = fullfile(dataPath, sub_str);
outPath = fullfile(basePath, 'LocalProcessed');
outSubPath = fullfile(outPath, sub_str);
if ~exist(outSubPath, 'dir')
    mkdir(outSubPath);
end

switch sub_str
    case 'S352'
        SETTING = 1;
        sess_files{1} = 'S352-022_SPECIALevents.mat';
    case 'S357'
        SETTING = 1;
        sess_files{1} = '357-027_SPECIALevents.mat';
        sess_files{2} = '357-028_SPECIALevents.mat';
        sess_files{3} = '357-029_SPECIALevents.mat';
    case 'S362'
        SETTING = 1;
        sess_files{1} = 'S362-055_SPECIALevents.mat';
        sess_files{2} = 'S362-056_SPECIALevents.mat';
        sess_files{3} = 'S362-057_SPECIALevents.mat';
        sess_files{4} = 'S362-058_SPECIALevents.mat';
    case 'S369'
        SETTING = 2;
        sess_files{1} = 'S369-043.mat';
        sess_files{2} = 'S369-044.mat';
        sess_files{3} = 'S369-046.mat';
        sess_files{4} = 'S369-049.mat';
    case 'S372'
        SETTING = 3;
        sess_files{1} = 'S372-024.mat';
        sess_files{2} = 'S372-025.mat';
        sess_files{3} = 'S372-026.mat';
        sess_files{4} = 'S372-027.mat';
    case 'S376'
        SETTING = 4;
        sess_files{1} = '376-039_Repetition Suppression.mat';
        sess_files{2} = '376-040_Repetition Suppression.mat';
        sess_files{3} = '376-041_Repetition Suppression.mat';
        sess_files{4} = '376-042_Repetition Suppression.mat';
end

ECoG_Fs = 2000;
Audio_Fs = 16000;
ECoGALL = Create_ECoG_Dataset(fullfile(subPath, sess_files), SETTING, ECoG_Fs, Audio_Fs);

for i=1:length(ECoGALL)
    ECoGALL{i} = Early_Preprocessing(ECoGALL{i}, line_noise_exit_flag);
end

if line_noise_exit_flag
    outFile = fullfile(outSubPath, 'ECoGALL_LinenoiseOnly.mat');
else
    outFile = fullfile(outSubPath, 'ECoGALL.mat');
end
save(outFile, 'ECoGALL', '-v7.3');

