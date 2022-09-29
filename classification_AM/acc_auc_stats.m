%%%% analyze accuracy and area-under-the-curve data from machine learning analysis
% data compiled by function /usr2/postdoc/amsmeier/compile_auc_and_acc.m
%
% updated by Andrew Meier 2022/2/2

auc_acc_data_file = '/usr2/postdoc/amsmeier/electrodes_acc_auc.mat';
 
close all
%  load(auc_acc_data_file); 
% rank_acc_auc(); % get perentile performance rankings
 
align_to_plot = 1; % plot onset-aligned
%     align_to_plot = 2; % plot stim-aligned

labeltype_to_plot = 'consonants';
%     labeltype_to_plot = 'vowel';
%     labeltype_to_plot = 'word';

metric_to_plot = 'acc'; % accuracy
% metric_to_plot = 'auc'; % area under the curve

clust_rank_thresh = 0.7; % ranges 0-1 (negative to not use); plot electrodes with this percentile performance or better

show_plots = 'on';
% show_plots = 'off';

 stim_key_file = '/projectnb2/busplab/Experiments/ECoG_Preprocessed_AM/stim_index_key.xlsx';
label_types = {'consonants', 'vowel', 'word'};
align_conds = {'onset', 'stim'}; 
sub_nums = [357; 362; 369; 372; 376]; 

nsubs = length(sub_nums);
n_align_conds = length(align_conds); 
n_label_types = length(label_types);
sub_names = {};
for isub = 1:nsubs
    sub_names{isub} = ['S', num2str(sub_nums(isub))];
end
nancol = NaN(height(electrodes), 1); 
stimkey = readtable(stim_key_file);

var_to_plot = [metric_to_plot, '_', labeltype_to_plot];

elc_to_plot = electrodes.type == align_to_plot; 
elc_to_plot = elc_to_plot & electrodes.(['clustrank_', var_to_plot]) > clust_rank_thresh; 

[p,~, stats] = anova1(electrodes.(var_to_plot)(elc_to_plot), electrodes.global_clust_num(elc_to_plot), show_plots)


