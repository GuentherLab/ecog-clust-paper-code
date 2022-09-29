% function multinode_feat_importance(temp_struct)

% params = load(temp_struct.temp_params_filename).params; 
% node_number = temp_struct.node_number; 

params = p; 
node_number = 5; 

% edb = load(fullfile(params.times_path, 'electrodes_database.mat')).electrodes_database;
node_electrode_rows = edb(edb.node_num == node_number, :); 

for electrode_idx = 1%:height(node_electrode_rows)
    sub_num = node_electrode_rows.subject(electrode_idx); 
    electrode_type = node_electrode_rows.type(electrode_idx); 
    electrode_id = node_electrode_rows.electrode(electrode_idx); 
    
    switch electrode_type
        case electrode_type == 2
            alignment = 'stim'; 
        case electrode_type == 1
            alignment = 'onset';
    end
    
    relevant_table = edb(edb.subject == sub_num & edb.type == electrode_type, :); 
    electrode_excluded_table = relevant_table(~(relevant_table.subject == sub_num & relevant_table.type == electrode_type & relevant_table.electrode == electrode_id), :);
    
    excluded_feat_struct = format_grouped_features(params, electrode_excluded_table);
    
    group_vals = fieldnames(excluded_feat_struct);
    
    for group_value_idx = 1:numel(group_vals)

    params.current_group_value = group_vals{group_value_idx};
    params.pooled_electrode_feature_set = excluded_feat_struct.(params.current_group_value).all_electrodes;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%   Begin the Training/Analysis   %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    excluded_electrode_checkpoint_file = sprintf('s%d_%s_elect_excl_%d_checkpoint_e%d_w%d_s%d.mat', sub_num, alignment, electrode_id, params.event_duration, params.window, params.stride);
    excluded_electrode_checkpoint_filename_full = fullfile(params.grouping_path, params.current_group_value, excluded_electrode_checkpoint_file);

    if ~isfile(excluded_electrode_checkpoint_filename_full)

        excluded_electrode_all_class_sub_cv_mRMR_LDA_data = struct();

        class_labels = params.class_names;
        for class_label_idx = 1:length(class_labels)
            class_label = class_labels{class_label_idx};
            fprintf('Class label: %s\n', class_label);

            excluded_electrode_all_class_sub_cv_mRMR_LDA_data.(class_label) = sub_spec_excluded_electrode_cv_mRMR_LDA(params, class_label); 

        end

        %%%% This is a useful block to have, next time you need it it won't take 15
        %%%% years to run. Just load final_checkpoint.mat
        fprintf('Saving checkpoint file... ');
        save(excluded_electrode_checkpoint_filename_full,...
            'excluded_electrode_all_class_sub_cv_mRMR_LDA_data',...
            '-v7.3');
        fprintf('Done.\n');

    else
        fprintf('Steps 1-4 already saved for %s.  Loading checkpoint file... ', params.current_group_value);
        load(excluded_electrode_checkpoint_filename_full,...
            'excluded_electrode_all_class_sub_cv_mRMR_LDA_data');

        fprintf('Done.\n');
    end
    
    end
    
    
end

% end 

