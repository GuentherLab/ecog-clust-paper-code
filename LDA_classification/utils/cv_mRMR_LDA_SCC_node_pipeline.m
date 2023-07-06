function cv_mRMR_LDA_SCC_node_pipeline(temp_struct)
% cv_mRMR_LDA_SCC_node_pipeline accepts the temporary struct defined in 
%   execute_multinode_pipeline. This function is a condensed version of the 
%   ML pipeline.
%
%   This version contains the processes that are necessary for an iterated
%   permutation test to randomize models for AUC statistics. The results
%   from 12 iterations (permutations) are saved in the auc_results_path.
%   The filenames contain the node number, which is used to track the
%   progress of the distributed processes. 
%   
%   See also: execute_multinode_pipeline,
%   indiv_electrode_cv_mRMR_LDA_OvR_rand 

params = load(temp_struct.temp_params_filename).params; 
node_number = temp_struct.node_number; 
% Function Pipeline for submission to multiple nodes on SCC frame. 
%
% Change relevant filepathes before use

% cd '/usr3/graduate/lcj126/Shortcut_scripts_LJ/CV_mRMR_LDA_randomization_optimization/';

tStart = tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%    File Path to Data    %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%    Define Parameters    %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class_labels = params.class_names;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%  Create/Load Feature Set  %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Building Cross-Validated, mRMR, LDA models trained on ECoG data. \nMultinodal Processing.\n\n');

sub_stim_e_id_auc_all_labels_all_iter = table();
sub_onset_e_id_auc_all_labels_all_iter = table();

nrepeat_min = 1;
nrepeat_max = 12;

for nrepeat = nrepeat_min:nrepeat_max   % iterate 12 times on each computer/node

    fprintf('%%%%%%%%%%%%\t  Node %d, Iteration %d  %%%%%%%%%%%%\n\n', node_number, nrepeat);
    
    rng(100 * node_number + nrepeat);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%   Begin the Training/Analysis   %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    all_class_sub_indiv_elect_cv_mRMR_LDA_OvR_data = struct();

    for class_label_idx = 1:length(class_labels)
        class_label = class_labels{class_label_idx};
%         params_struct.class_label = class_label;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% indiv_elect_cv_mRMR_LDA_OvR.m
        fprintf('Data from individual electrodes, all subjects. One vs Rest\n\n');

        sub_indiv_elect_cv_mRMR_LDA_OvR_data_rand = indiv_electrode_cv_mRMR_LDA_OvR_rand(params, class_label);

        all_class_sub_indiv_elect_cv_mRMR_LDA_OvR_data.(class_label) = sub_indiv_elect_cv_mRMR_LDA_OvR_data_rand;
    end

    fprintf('CV, mRMR, LDA models from individual electrodes, all subjects successful.\n\nCalculating mean AUCs... ');
    ind_ovr = all_class_sub_indiv_elect_cv_mRMR_LDA_OvR_data;

    sub_nums = [357, 362, 369, 372, 376];

    sub_stim_e_id_auc_all_labels = table();

    for sub_num = sub_nums
        sub_str = strcat('S', string(sub_num));

        sub_stim_e_ids = ind_ovr.vowel_name.(sub_str).stim_data.stim_e_id; 

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

        sub_onset_e_ids = ind_ovr.vowel_name.(sub_str).onset_data.onset_e_id; 

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

    stim_cluster_col = cell([height(sub_stim_e_id_auc_all_labels), 1]);
    stim_roi_col = cell([height(sub_stim_e_id_auc_all_labels), 1]);

    for i = 1:height(sub_stim_e_id_auc_all_labels)
        sub_num = sub_stim_e_id_auc_all_labels.('sub_num')(i);
        e_id = sub_stim_e_id_auc_all_labels.('stim_e_id')(i);

        e_cluster = get_cluster_name(sub_num, e_id, 'stim');
        e_roi = get_ROI_name(sub_num, e_id, 'stim');

        stim_cluster_col{i} = e_cluster;
        stim_roi_col{i} = e_roi;
    end
    stim_cluster_col = categorical(stim_cluster_col);
    stim_roi_col = categorical(stim_roi_col);
    sub_stim_e_id_auc_all_labels = addvars(sub_stim_e_id_auc_all_labels,...
        stim_cluster_col, stim_roi_col,...
        'After', 'stim_e_id',...
        'NewVariableNames', {'stim_cluster', 'stim_roi'});

    sub_stim_e_id_auc_all_labels_all_iter = [sub_stim_e_id_auc_all_labels_all_iter; sub_stim_e_id_auc_all_labels]; 
    
    onset_cluster_col = cell([height(sub_onset_e_id_auc_all_labels), 1]);
    onset_roi_col = cell([height(sub_onset_e_id_auc_all_labels), 1]);

    for i = 1:height(sub_onset_e_id_auc_all_labels)
        sub_num = sub_onset_e_id_auc_all_labels.('sub_num')(i);
        e_id = sub_onset_e_id_auc_all_labels.('onset_e_id')(i);

        e_cluster = get_cluster_name(sub_num, e_id, 'onset');
        e_roi = get_ROI_name(sub_num, e_id, 'onset');

        onset_cluster_col{i} = e_cluster;
        onset_roi_col{i} = e_roi;        
    end
    onset_cluster_col = categorical(onset_cluster_col);
    onset_roi_col = categorical(onset_roi_col);
    sub_onset_e_id_auc_all_labels = addvars(sub_onset_e_id_auc_all_labels,...
        onset_cluster_col, onset_roi_col,...
        'After', 'onset_e_id',...
        'NewVariableNames', {'onset_cluster', 'onset_roi'});

    sub_onset_e_id_auc_all_labels_all_iter = [sub_onset_e_id_auc_all_labels_all_iter; sub_onset_e_id_auc_all_labels]; 
    
    fprintf('Done.\n');
end
tEnd = toc(tStart);
fprintf('All iterations, node %d\n\nElapsed time: %f. \n\nSaving...', node_number, tEnd);
mkdir(fullfile(params.grouping_path, params.current_group_value,'AUC_Results'));

all_nodes_auc_res_file = sprintf('node_%d_iter_%d_thru_%d_indiv_elect_mean_auc', node_number, nrepeat_min, nrepeat_max);
all_nodes_auc_res_filename_full = fullfile(params.grouping_path, params.current_group_value,'AUC_Results', all_nodes_auc_res_file);
save(all_nodes_auc_res_filename_full, 'sub_stim_e_id_auc_all_labels_all_iter', 'sub_onset_e_id_auc_all_labels_all_iter', '-v7.3');

fprintf('Done.\n\n');

end














