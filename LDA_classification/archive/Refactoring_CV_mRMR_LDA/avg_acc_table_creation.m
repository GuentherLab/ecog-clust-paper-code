class_labels = fieldnames(all_class_sub_cv_mRMR_LDA_data); 
avg_acc_table = table();
for class_label_idx = 1:length(class_labels)
    class_label = class_labels{class_label_idx} 
    class_data = all_class_sub_cv_mRMR_LDA_data.(class_label);
    
    sub_strs = fieldnames(class_data); 
    
    t_temp = table();
    for sub_str_idx = 1:length(sub_strs)
        sub_str = sub_strs{sub_str_idx}; 
        [s, e] = regexp(sub_str, '(\d+)');
        sub_num = str2double(sub_str(s:e)) 
        sub_class_data = class_data.(sub_str); 
        
        sub_class_stim_data = sub_class_data.stim_data; 
        sub_class_onset_data = sub_class_data.onset_data; 
        
        mean_stim_acc = mean(sub_class_stim_data.(2).stim_acc)
        mean_onset_acc = mean(sub_class_onset_data.(2).onset_acc)
        
        % choosing to do onset in first column and stim second because of
        % type number in edb (stim = 2, onset = 1)
        t_temp = [t_temp; table(table(mean_onset_acc, mean_stim_acc,...
            'VariableNames', {'mean_onset_acc', 'mean_stim_acc'}),...
            'VariableNames', {class_label})];
        
    end
    avg_acc_table = [avg_acc_table, t_temp];
end
avg_acc_table = [table(p.sub_nums', 'VariableNames', {'sub_num'}), avg_acc_table]
save(fullfile(p.grouping_path, p.current_group_value, 'avg_acc_table_e2000_w50_s25.mat'), 'avg_acc_table');

%%
load(fullfile(p.grouping_path, p.current_group_value, 'avg_acc_table_e2000_w50_s25.mat'), 'avg_acc_table');


