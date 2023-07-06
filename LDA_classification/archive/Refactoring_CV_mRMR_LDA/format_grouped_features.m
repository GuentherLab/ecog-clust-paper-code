function feat_struct = format_grouped_features(params, edb)
% format_grouped_features Formats a struct based on the grouping variable
%   in parametersClass object. 
%   Output struct has fields corresponding to unique values of the grouping variable. 

sub_nums = unique(edb.subject); 
class_names = params.class_names; 
grp_var = params.grouping_variable;

if ~strcmp(grp_var, 'none')
    unique_grp_vals = unique(edb.(grp_var));
else
    unique_grp_vals = {'none'};
end

feat_struct = struct();
    for grp_val_idx = 1:numel(unique_grp_vals)
        all_subs_feat_table = table(sub_nums, 'VariableNames', {'sub_num'});
        all_subs_indiv_elect_feat_table = table(sub_nums, 'VariableNames', {'sub_num'});

        grp_val = unique_grp_vals{grp_val_idx}; 

        if strcmp(grp_val, 'none')
            relevant_electrodes_table = edb;
        else
            relevant_electrodes_table = edb(strcmp(edb.(grp_var), grp_val), :);
        end

        for class_name_idx = 1:numel(class_names)
            class_name = class_names{class_name_idx};
            indiv_electrodes_class_name_table = table();
            pooled_electrodes_class_name_table = table();
            for sub_num_idx = 1:numel(sub_nums)
                sub_num = sub_nums(sub_num_idx);

                class_labels_obj = classLabelsClass(sub_num);
                class_labels = class_labels_obj.(class_name);
                class_labels_ohe = class_labels_obj.([class_name, '_ohe']);

                relevant_sub_grp_electrodes = relevant_electrodes_table(relevant_electrodes_table.subject == sub_num, :);

                stim_pooled_electrodes = table();
                onset_pooled_electrodes = table();
                for relevant_electrodes_row = 1:height(relevant_sub_grp_electrodes)
                    feat_table_temp = relevant_sub_grp_electrodes.feat_set{relevant_electrodes_row}; 
                    if relevant_sub_grp_electrodes.type(relevant_electrodes_row) == 2
                        relevant_sub_grp_electrodes.feat_set{relevant_electrodes_row} = table(feat_table_temp, class_labels, class_labels_ohe,...
                            'VariableNames', {'stim_feat', 'class_labels', 'class_labels_ohe'});
                        stim_pooled_electrodes = [stim_pooled_electrodes, feat_table_temp];
                    else
                        relevant_sub_grp_electrodes.feat_set{relevant_electrodes_row} = table(feat_table_temp, class_labels, class_labels_ohe,...
                            'VariableNames', {'onset_feat', 'class_labels', 'class_labels_ohe'});
                        onset_pooled_electrodes = [onset_pooled_electrodes, feat_table_temp];
                    end
                end

                if sum(size(stim_pooled_electrodes)) == 0
                    stim_pooled_electrodes = table(zeros([height(class_labels), 1])); 
                elseif sum(size(onset_pooled_electrodes)) == 0
                    onset_pooled_electrodes = table(zeros([height(class_labels), 1])); 
                end

                stim_pooled_electrodes = table(stim_pooled_electrodes, class_labels, class_labels_ohe, 'VariableNames', {'stim_feat', 'class_labels', 'class_labels_ohe'});
                onset_pooled_electrodes = table(onset_pooled_electrodes, class_labels, class_labels_ohe, 'VariableNames', {'onset_feat', 'class_labels', 'class_labels_ohe'});

                pooled_electrodes_class_name_table_temp = table({stim_pooled_electrodes}, {onset_pooled_electrodes},...
                    'VariableNames', {'stim_feat', 'onset_feat'});
                pooled_electrodes_class_name_table = [pooled_electrodes_class_name_table; pooled_electrodes_class_name_table_temp];

                stim_sub_grp_electrodes = relevant_sub_grp_electrodes(relevant_sub_grp_electrodes.type == 2,{'electrode', 'feat_set'});
                stim_sub_grp_electrodes = renamevars(stim_sub_grp_electrodes, {'electrode', 'feat_set'}, {'stim_e_id', 'stim_feat_class'});
                onset_sub_grp_electrodes = relevant_sub_grp_electrodes(relevant_sub_grp_electrodes.type == 1,{'electrode', 'feat_set'});
                onset_sub_grp_electrodes = renamevars(onset_sub_grp_electrodes, {'electrode', 'feat_set'}, {'onset_e_id', 'onset_feat_class'});

                indiv_electrodes_class_name_table_temp = table({stim_sub_grp_electrodes}, {onset_sub_grp_electrodes},...
                    'VariableNames', {'stim_feat', 'onset_feat'});
                indiv_electrodes_class_name_table = [indiv_electrodes_class_name_table; indiv_electrodes_class_name_table_temp];

            end

            all_subs_feat_table = [all_subs_feat_table, table(pooled_electrodes_class_name_table, 'VariableNames', {class_name})];
            all_subs_indiv_elect_feat_table = [all_subs_indiv_elect_feat_table, table(indiv_electrodes_class_name_table, 'VariableNames', {class_name})];

        end
        feat_struct.(grp_val).all_electrodes = all_subs_feat_table;
        feat_struct.(grp_val).indiv_electrodes = all_subs_indiv_elect_feat_table;

    end
end

