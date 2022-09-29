clear all; close all; clc;

addpath('PreprocessTools');
addpath(genpath('EEG-Clean-Tools-master'));
addpath(genpath('eeglab13_5_4b'));
    
sub_str = 'S376'; % {'S352', 'S357', 'S362', 'S369', 'S372', 'S376'}
filt_type = {'Hilbert'}; % {'Hilbert, 'MultiTaper', 'Filter', 'Wavelet'};
processing_type = 'highgamma';

switch processing_type
    case 'highgamma'
        save_label = 'HG';
    case 'beta'
        save_label = 'BT';
end


basePath = error('Enter data directory here');
outPath = fullfile(basePath, 'LocalProcessed');
outSubPath = fullfile(outPath, sub_str);

inFile = fullfile(outSubPath, 'ECoGALL.mat');
load(inFile);

if any(cellfun(@(x) strcmp(x, 'Hilbert'), filt_type))

    for i=1:length(ECoGALL)
        ECoGALL{i} = ECoG_Hilbert_HG(ECoGALL{i}, processing_type);
    end

    outFile = fullfile(outSubPath, ['ECoGALL_Hilbert_' save_label '.mat']);
    save(outFile, 'ECoGALL', '-v7.3');

end

if any(cellfun(@(x) strcmp(x, 'MultiTaper'), filt_type))

    for i=1:length(ECoGALL)
        ECoGALL{i} = ECoG_Multi_Taper(ECoGALL{i}, processing_type);
    end

    outFile = fullfile(outSubPath, ['ECoGALL_MultiTaper_' save_label '.mat']);
    save(outFile, 'ECoGALL', '-v7.3');

end

if any(cellfun(@(x) strcmp(x, 'Filter'), filt_type))
    
   for i=1:length(ECoGALL)
        ECoGALL{i} = ECoG_Filter(ECoGALL{i}, processing_type);
    end

    outFile = fullfile(outSubPath, ['ECoGALL_Filter_' save_label '.mat']);
    save(outFile, 'ECoGALL', '-v7.3');

end 

if any(cellfun(@(x) strcmp(x, 'Wavelet'), filt_type))
    
	for i=1:length(ECoGALL)
        ECoGALL{i} = ECoG_Wavelet(ECoGALL{i}, processing_type);
	end

    outFile = fullfile(outSubPath, ['ECoGALL_Wavelet_' save_label '.mat']);
    save(outFile, 'ECoGALL', '-v7.3');

end 