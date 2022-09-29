  %%% load electrodes table with leave-one-out coding scores for word, consonant, vowel
 % adapted from /usr2/postdoc/amsmeier/top_raw_pct_p_val_auto_AM.m
 %
 %%% updated 2022/3/9 by AM
 
 clear % if commented in, this script will always reload data, rather than using existing data in the workspace
 
 
 save_data = 0;
    savename = '/usr2/postdoc/amsmeier/leave_one_out_data';
 
 
 %% this section takes a few minutes to run
 if ~exist('p', 'var')
    p = parametersClass('/projectnb/busplab/Experiments/ECoG_fMRI_RS/Experiments/ECoG_Preprocessed_LJ/',... % path to data
        2000,...            % event duration (ms)
        50,...              % window (ms)
        25,...              % stride (ms)
        'none',...    % grouping variable
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
all_stim_counts = groupsummary(edb_nonan(edb_nonan.type == 2, :), {'cluster_name'},...
    'IncludeEmptyGroups', true, 'IncludeMissingGroups', true);
all_onset_counts = groupsummary(edb_nonan(edb_nonan.type == 1, :), {'cluster_name'},...
    'IncludeEmptyGroups', true, 'IncludeMissingGroups', true);

class_labels = p.class_names; 
data_cols = {'word_accuracy_change_wo_electrode',...
    'cons_accuracy_change_wo_electrode',...
    'vowel_accuracy_change_wo_electrode'};

%%
plotops.include_unused_clusters = 0; % include or don't include unused clusters in the bar plots
plotops.subjects_as_colors = 0; % separate out bars by subjects' contributions


per_sub_flag = false;
swc = table();
swc.cluster_name = unique(edb_nonan.cluster_name);
scc = swc;
svc = swc;

owc = swc; 
occ = swc; 
ovc = swc;

class_label = 'word_name';
[sw_sorted, ow_sorted] = sort_by_score(edb_nonan, class_label);


%
if save_data
    save(savename)
end




%%
function [sub_stim_sorted, sub_onset_sorted] = sort_by_score(edb_nonan, class_label)

data_cols = {'word_accuracy_change_wo_electrode',...
    'cons_accuracy_change_wo_electrode',...
    'vowel_accuracy_change_wo_electrode'};
switch class_label
    case 'word_name'
        data_col = 1; 
    case 'consonants_name'
        data_col = 2; 
    case 'vowel_name'
        data_col = 3; 
end

sub_stim_sorted = sortrows(edb_nonan(edb_nonan.type == 2,...
    :),...
    {data_cols{data_col}},...
    {'descend'});
sub_onset_sorted = sortrows(edb_nonan(edb_nonan.type == 1,...
    :),...
    {data_cols{data_col}},...
    {'descend'});

end