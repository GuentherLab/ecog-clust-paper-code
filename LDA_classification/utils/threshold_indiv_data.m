function threshd_auc_data = threshold_indiv_data(sub_stim_e_id_auc_all_labels, sub_onset_e_id_auc_all_labels, threshold_table)

ovr_labels = threshold_table.Properties.VariableNames;

threshd_auc_data = table(); 

for ovr_label_idx = 1:numel(ovr_labels)
    ovr_label = ovr_labels{ovr_label_idx}; 
    
    stim_ovr_label_auc = sub_stim_e_id_auc_all_labels.(ovr_label); 
    stim_sig_e_data = sub_stim_e_id_auc_all_labels(stim_ovr_label_auc >= threshold_table{'stim', ovr_label}, {'sub_num', 'stim_e_id', 'stim_table_idx', 'stim_cluster', 'stim_roi', ovr_label});
    stim_sig_e_data.Properties.VariableNames{ovr_label} = 'mean_auc'; 
    
    onset_ovr_label_auc = sub_onset_e_id_auc_all_labels.(ovr_label); 
    onset_sig_e_data = sub_onset_e_id_auc_all_labels(onset_ovr_label_auc >= threshold_table{'onset', ovr_label}, {'sub_num', 'onset_e_id', 'onset_table_idx', 'onset_cluster', 'onset_roi', ovr_label});
    onset_sig_e_data.Properties.VariableNames{ovr_label} = 'mean_auc';     
    
    threshd_auc_data_temp = table(categorical({ovr_label}), {stim_sig_e_data}, {onset_sig_e_data}, 'VariableNames', {'ovr_label', 'stim_sig_data', 'onset_sig_data'});
    threshd_auc_data = [threshd_auc_data; threshd_auc_data_temp];
    
end

% i = 1; j = 1; k = 1; 
% for ovr_label_idx = 1:numel(ovr_labels)
%     ovr_label = ovr_labels(ovr_label_idx); 
%     
%     if length(ovr_label) == 4
%         word_ovr{i} = ovr_label; 
%         i = i+1; 
%     elseif length(ovr_label) == 2
%         cons_ovr{j} = ovr_label;
%         j = j+1;
%     elseif length(ovr_label) == 1
%         vowel_ovr{k} = ovr_label;
%         k = k+1; 
%     end
% 
% end

end









