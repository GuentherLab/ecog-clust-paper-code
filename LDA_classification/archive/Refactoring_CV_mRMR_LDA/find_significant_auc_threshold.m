function threshold_table = find_significant_auc_threshold(stim_data_all_iterations, onset_data_all_iterations, p_value)
% Returns table of minimum AUC values (thresholds).
% 
% threshold_table = find_significant_auc_threshold(stim_data_all_iterations, onset_data_all_iterations, p_value)
% 
%   Returns a table of mean AUC values that serve as minimum 
%   thresholds for determining significant electrodes base on real data.
%   
%   Expects mean AUC tables (both alignment cases) resulting from 
%   execute_multinode_pipeline.m

z_thresh = icdf('normal', 1-p_value, 0, 1); 

ovr_labels = stim_data_all_iterations.Properties.VariableNames; 
ovr_labels(cellfun(@(x) strcmp(x, 'sub_num') || strcmp(x, 'stim_e_id') || strcmp(x, 'stim_cluster') || strcmp(x, 'stim_roi') , ovr_labels)) = [];

stim_pop_norm_min_auc_thresh = table();
onset_pop_norm_min_auc_thresh = table();

for ovr_label_idx = 1:length(ovr_labels)
    ovr_label = ovr_labels{ovr_label_idx};

%%%% STIM
    
    stim_sorted_mean_auc = sortrows(stim_data_all_iterations, ovr_label, 'ascend');
    stim_ovr_label_mean_auc = stim_sorted_mean_auc(:,{'sub_num', 'stim_e_id', 'stim_cluster', 'stim_roi', ovr_label});

    stim_probs = stim_ovr_label_mean_auc.(ovr_label); 
    stim_z_scores = zscore(stim_probs);
    
    stim_sig_e_data = stim_ovr_label_mean_auc(stim_z_scores >= z_thresh, :);
    stim_pop_norm_min_auc_thresh_temp = table(min(stim_sig_e_data.(ovr_label)), 'VariableNames', {ovr_label});
    stim_pop_norm_min_auc_thresh = [stim_pop_norm_min_auc_thresh, stim_pop_norm_min_auc_thresh_temp];
    
%%%% ONSET

    onset_sorted_mean_auc = sortrows(onset_data_all_iterations, ovr_label, 'ascend');
    onset_ovr_label_mean_auc = onset_sorted_mean_auc(:,{'sub_num', 'onset_e_id', 'onset_cluster', 'onset_roi', ovr_label});

    onset_probs = onset_ovr_label_mean_auc.(ovr_label); 
    onset_z_scores = zscore(onset_probs);

    onset_sig_e_data = onset_ovr_label_mean_auc(onset_z_scores >= z_thresh, :);
    onset_pop_norm_min_auc_thresh_temp = table(min(onset_sig_e_data.(ovr_label)), 'VariableNames', {ovr_label});
    onset_pop_norm_min_auc_thresh = [onset_pop_norm_min_auc_thresh, onset_pop_norm_min_auc_thresh_temp];

end

threshold_table = [stim_pop_norm_min_auc_thresh; onset_pop_norm_min_auc_thresh];
threshold_table.Properties.RowNames = {'stim', 'onset'};

end




