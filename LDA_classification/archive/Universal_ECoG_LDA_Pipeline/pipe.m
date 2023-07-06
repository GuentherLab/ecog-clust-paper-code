close all; 
clear all; 
clc;

tStart = tic;
p = parametersClass('/projectnb/busplab/Experiments/ECoG_fMRI_RS/Experiments/ECoG_Preprocessed_LJ/',... % path to data
    2000,...            % event duration (ms)
    50,...              % window (ms)
    25,...              % stride (ms)
    'none',...    % grouping variable
    150,...             % topN feat pooled
    30)                 % topN feat indiv

[db, edb] = p.generate_database; 

grouped_feat_struct = format_grouped_features(p, edb); 

group_vals = fieldnames(grouped_feat_struct);


