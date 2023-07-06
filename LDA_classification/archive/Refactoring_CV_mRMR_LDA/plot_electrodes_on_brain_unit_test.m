% run this section to generate/load the data 

clc
set(0,'DefaultFigureWindowStyle','docked')
addpath /project/busplab/software/ecog/scripts_LJ/Refactoring_CV_mRMR_LDA

if ~exist('p', 'var')
    p = parametersClass('/projectnb/busplab/Experiments/ECoG_fMRI_RS/Experiments/ECoG_Preprocessed_LJ/',... % path to data
        2000,...            % event duration (ms)
        50,...              % window (ms)
        25,...              % stride (ms)
        'none',...          % grouping variable
        150,...             % topN feat pooled
        30)                 % topN feat indiv
end
if ~exist('edb', 'var') 
    [~, edb] = p.generate_database; 
end
if ~exist('edb_nonan', 'var') 
    generate_edb_nonan
end

%%
close all force;
% change this subset to control which electrodes get plotted:
% 2 = stim, 1 = onset
alignment = 1; 

% {'ESP-s'} = 1
% {'ESP-r'} = 2
% {'PtM-s'} = 3
% {'PtM-r'} = 4
% {'ME-sb'} = 5
% {'ME-sn'} = 6
% {'AP-r' } = 7
% {'AP-s' } = 8
global_cluster_number = 6; % clusterkey handles labeling in plot fxn

edb_nonan_subset = edb_nonan(edb_nonan.type == alignment & edb_nonan.global_clust_num == global_cluster_number,:);

% some options
ops = struct(); 
ops.force_close = 1; 
ops.clusters_as_colors = 1; 
ops.subjects_as_shapes = 1; 
ops.suppress_brain_surfs = 0; 
ops.do_mosaic = 0;
ops.selected_hemi = 'left';
ops.viewpoint = 'left';
% ops.img_savename = 'test.png';

t = plot_electrodes_on_brain(edb_nonan_subset, 'word_name', ops);

% this is functioning as my shape legend for now:
unique(t)
