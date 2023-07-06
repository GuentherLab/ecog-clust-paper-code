
function multinode_feat_importance(temp_struct)

params = load(temp_struct.temp_params_filename).params; 
node_number = temp_struct.node_number; 

% params = p; 
% node_number = 26; 

edb = load(fullfile(params.times_path, 'electrodes_database.mat')).electrodes_database;
node_electrode_rows = edb(edb.node_num == node_number, :); 

for electrode_idx = 1:height(node_electrode_rows)
    sub_num = node_electrode_rows.subject(electrode_idx); 
    electrode_type = node_electrode_rows.type(electrode_idx);
    electrode_id = node_electrode_rows.electrode(electrode_idx);
    edb_row_num = node_electrode_rows.edb_row_num(electrode_idx);
    
    switch electrode_type
        case 2
            alignment = 'stim'; 
        case 1
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

        excluded_electrode_checkpoint_file = sprintf('row%d_s%d_%s_elect_excl_%d_checkpoint_e%d_w%d_s%d.mat', edb_row_num, sub_num, alignment, electrode_id, params.event_duration, params.window, params.stride);
        excluded_electrode_checkpoint_filename_full = fullfile(params.grouping_path, params.current_group_value, 'excluded_electrode_data', excluded_electrode_checkpoint_file);

        if ~isfile(excluded_electrode_checkpoint_filename_full)

            excluded_electrode_all_class_sub_cv_mRMR_LDA_data = struct();

            class_labels = params.class_names;
            for class_label_idx = 1:length(class_labels)
                class_label = class_labels{class_label_idx};

                excluded_electrode_all_class_sub_cv_mRMR_LDA_data.(class_label) = test_LDA_fxn(params, class_label);

            end

%             disp(excluded_electrode_all_class_sub_cv_mRMR_LDA_data)

            save(excluded_electrode_checkpoint_filename_full, '-struct', 'excluded_electrode_all_class_sub_cv_mRMR_LDA_data', '-v7.3'); 

%         else
%             load(excluded_electrode_checkpoint_filename_full, 'excluded_electrode_all_class_sub_cv_mRMR_LDA_data');
        end

%         edb_matfile = matfile(fullfile(params.times_path, 'electrodes_database.mat'), 'Writable', true); 
%         electrodes_database = edb_matfile.electrodes_database; 
%         electrodes_database.word_accuracy_wo_electrode((electrodes_database.subject == sub_num &...
%             electrodes_database.type == electrode_type &...
%             electrodes_database.electrode == electrode_id)) = excluded_electrode_all_class_sub_cv_mRMR_LDA_data.word_name; 
% 
%         electrodes_database.cons_accuracy_wo_electrode((electrodes_database.subject == sub_num &...
%             electrodes_database.type == electrode_type &...
%             electrodes_database.electrode == electrode_id)) = excluded_electrode_all_class_sub_cv_mRMR_LDA_data.consonants_name; 
% 
%         electrodes_database.vowel_accuracy_wo_electrode((electrodes_database.subject == sub_num &...
%             electrodes_database.type == electrode_type &...
%             electrodes_database.electrode == electrode_id)) = excluded_electrode_all_class_sub_cv_mRMR_LDA_data.vowel_name; 
% 
%         edb_matfile.electrodes_database = electrodes_database;

    end
end

end 

