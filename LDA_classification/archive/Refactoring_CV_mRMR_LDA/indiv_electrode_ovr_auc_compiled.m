% indiv_electrode_ovr_auc_compiled operates on individual electrode AUC
% results. Note depending on grouping variable chose, you may have to add
% the variable of interest to the bottom where results are compiled into
% the output table. (lines 111 and 138)

ind_ovr = all_class_sub_indiv_elect_cv_mRMR_LDA_OvR_data;

class_labels = fieldnames(ind_ovr); 

sub_nums = [357, 362, 369, 372, 376];

sub_stim_e_id_auc_all_labels = table();

for sub_num = sub_nums
    sub_str = strcat('S', string(sub_num));
    
    sub_stim_e_ids = ind_ovr.word_name.(sub_str).stim_data.stim_e_id; 

    for stim_e_id_idx = 1:length(sub_stim_e_ids)
        sub_stim_e_id = sub_stim_e_ids(stim_e_id_idx);

        stim_ovr_label_mean_auc_temp_temp = table();
        
        for class_label_idx = 1:length(class_labels)
            class_label = class_labels{class_label_idx};
            
            sub_stim_class_e_id_data = ind_ovr.(class_label).(sub_str).stim_data.stim_fold_label_table{stim_e_id_idx};
            
            ovr_labels = sub_stim_class_e_id_data.Properties.VariableNames; 
            ovr_labels(cellfun(@(x) strcmp(x, 'fold'), ovr_labels)) = [];

            for ovr_label_idx = 1:length(ovr_labels)
                ovr_label = ovr_labels{ovr_label_idx}; 
                
                stim_ovr_label_mean_auc = mean(sub_stim_class_e_id_data.(ovr_label).stim_auc);
                stim_ovr_label_mean_auc_temp = table(stim_ovr_label_mean_auc, 'VariableNames', {ovr_label});
                
                stim_ovr_label_mean_auc_temp_temp = [stim_ovr_label_mean_auc_temp_temp, stim_ovr_label_mean_auc_temp];
                
            end
            stim_ovr_label_mean_auc_temp_temp_temp = [table(sub_num, 'VariableNames', {'sub_num'}), table(sub_stim_e_id, 'VariableNames', {'stim_e_id'}), stim_ovr_label_mean_auc_temp_temp];
            
        end
        sub_stim_e_id_auc_all_labels = [sub_stim_e_id_auc_all_labels; stim_ovr_label_mean_auc_temp_temp_temp];
        
    end
end

sub_onset_e_id_auc_all_labels = table();

for sub_num = sub_nums
    sub_str = strcat('S', string(sub_num));
    
    sub_onset_e_ids = ind_ovr.word_name.(sub_str).onset_data.onset_e_id; 

    for onset_e_id_idx = 1:length(sub_onset_e_ids)
        sub_onset_e_id = sub_onset_e_ids(onset_e_id_idx);

        onset_ovr_label_mean_auc_temp_temp = table();
        
        for class_label_idx = 1:length(class_labels)
            class_label = class_labels{class_label_idx};
            
            sub_onset_class_e_id_data = ind_ovr.(class_label).(sub_str).onset_data.onset_fold_label_table{onset_e_id_idx};
            
            ovr_labels = sub_onset_class_e_id_data.Properties.VariableNames; 
            ovr_labels(cellfun(@(x) strcmp(x, 'fold'), ovr_labels)) = [];

            for ovr_label_idx = 1:length(ovr_labels)
                ovr_label = ovr_labels{ovr_label_idx}; 
                
                onset_ovr_label_mean_auc = mean(sub_onset_class_e_id_data.(ovr_label).onset_auc);
                onset_ovr_label_mean_auc_temp = table(onset_ovr_label_mean_auc, 'VariableNames', {ovr_label});
                
                onset_ovr_label_mean_auc_temp_temp = [onset_ovr_label_mean_auc_temp_temp, onset_ovr_label_mean_auc_temp];
                
            end
            onset_ovr_label_mean_auc_temp_temp_temp = [table(sub_num, 'VariableNames', {'sub_num'}), table(sub_onset_e_id, 'VariableNames', {'onset_e_id'}), onset_ovr_label_mean_auc_temp_temp];
            
        end
        sub_onset_e_id_auc_all_labels = [sub_onset_e_id_auc_all_labels; onset_ovr_label_mean_auc_temp_temp_temp];
        
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%

stim_table_row_col = zeros([height(sub_stim_e_id_auc_all_labels), 1]);
stim_cluster_col = cell([height(sub_stim_e_id_auc_all_labels), 1]);
stim_roi_col = cell([height(sub_stim_e_id_auc_all_labels), 1]);
stim_depth_col = cell([height(sub_stim_e_id_auc_all_labels), 1]);

for i = 1:height(sub_stim_e_id_auc_all_labels)
    sub_num = sub_stim_e_id_auc_all_labels.('sub_num')(i);
    e_id = sub_stim_e_id_auc_all_labels.('stim_e_id')(i);
    
    e_table_row = get_electrode_table_row_idx(sub_num, e_id, 'stim');
    e_cluster = get_cluster_name(sub_num, e_id, 'stim');
    e_roi = get_ROI_name(sub_num, e_id, 'stim');
    e_depth = get_depth(sub_num, e_id, 'stim');
    
    stim_table_row_col(i) = e_table_row; 
    stim_cluster_col{i} = e_cluster;
    stim_roi_col{i} = e_roi;
    stim_depth_col{i} = e_depth;
end
stim_cluster_col = categorical(stim_cluster_col);
stim_roi_col = categorical(stim_roi_col);
stim_depth_col = categorical(stim_depth_col);
sub_stim_e_id_auc_all_labels = addvars(sub_stim_e_id_auc_all_labels,...
    stim_table_row_col, stim_cluster_col, stim_roi_col, stim_depth_col,...
    'After', 'stim_e_id',...
    'NewVariableNames', {'stim_table_idx', 'stim_cluster', 'stim_roi', 'stim_depth'});

onset_table_row_col = zeros([height(sub_onset_e_id_auc_all_labels), 1]);
onset_cluster_col = cell([height(sub_onset_e_id_auc_all_labels), 1]);
onset_roi_col = cell([height(sub_onset_e_id_auc_all_labels), 1]);
onset_depth_col = cell([height(sub_onset_e_id_auc_all_labels), 1]);

for i = 1:height(sub_onset_e_id_auc_all_labels)
    sub_num = sub_onset_e_id_auc_all_labels.('sub_num')(i);
    e_id = sub_onset_e_id_auc_all_labels.('onset_e_id')(i);
    
    e_table_row = get_electrode_table_row_idx(sub_num, e_id, 'onset');
    e_cluster = get_cluster_name(sub_num, e_id, 'onset');
    e_roi = get_ROI_name(sub_num, e_id, 'onset');
    e_depth = get_depth(sub_num, e_id, 'onset');
    
    onset_table_row_col(i) = e_table_row; 
    onset_cluster_col{i} = e_cluster;
    onset_roi_col{i} = e_roi;
    onset_depth_col{i} = e_depth;
end
onset_cluster_col = categorical(onset_cluster_col);
onset_roi_col = categorical(onset_roi_col);
onset_depth_col = categorical(onset_depth_col);
sub_onset_e_id_auc_all_labels = addvars(sub_onset_e_id_auc_all_labels,...
    onset_table_row_col, onset_cluster_col, onset_roi_col, onset_depth_col,...
    'After', 'onset_e_id',...
    'NewVariableNames', {'onset_table_idx', 'onset_cluster', 'onset_roi', 'onset_depth'});










