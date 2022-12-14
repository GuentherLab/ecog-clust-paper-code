close all; 
clear all; 
clc;

tStart = tic;
p = parametersClass('/projectnb/busplab/Experiments/ECoG_fMRI_RS/Experiments/ECoG_Preprocessed_LJ/',... % path to data
    2000,...            % event duration (ms)
    50,...              % window (ms)
    25,...              % stride (ms)
    'none',...    % grouping variable
    150,...             % topN feat pooled
    30)                 % topN feat indiv

[db, edb] = p.generate_database; 

grouped_feat_struct = format_grouped_features(p, edb); 

group_vals = fieldnames(grouped_feat_struct);
%%
for group_value_idx = 1:numel(group_vals)

    p.current_group_value = group_vals{group_value_idx};
    p.pooled_electrode_feature_set = grouped_feat_struct.(p.current_group_value).all_electrodes;
    p.individual_electrode_feature_set = grouped_feat_struct.(p.current_group_value).indiv_electrodes;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%   Begin the Training/Analysis   %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    checkpoint_file = sprintf('checkpoint_e%d_w%d_s%d.mat', p.event_duration, p.window, p.stride);
    checkpoint_filename_full = fullfile(p.grouping_path, p.current_group_value, checkpoint_file);

    if ~isfile(checkpoint_filename_full)

        all_class_sub_cv_mRMR_LDA_OvR_data = struct();
        all_class_sub_cv_mRMR_LDA_data = struct();
        all_class_sub_indiv_elect_cv_mRMR_LDA_OvR_data = struct();
        all_class_sub_indiv_elect_cv_mRMR_LDA_data = struct();

        class_labels = p.class_names;
        for class_label_idx = 1:length(class_labels)
            class_label = class_labels{class_label_idx};
            fprintf('Class label: %s\n', class_label);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            %%%% Step 1. all_electrodes_cv_mRMR_LDA_OvR.m
                fprintf('Step 1. All Electrodes, One vs. Rest\n');
            all_class_sub_cv_mRMR_LDA_OvR_data.(class_label) = all_electrodes_cv_mRMR_LDA_OvR(p, class_label);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% Step 2. all_elect_cv_mRMR_LDA.m
                fprintf('Step 2. All Electrodes, Regular Labels\n');
            all_class_sub_cv_mRMR_LDA_data.(class_label) = all_electrodes_cv_mRMR_LDA(p, class_label); 

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% Step 3. indiv_elect_cv_mRMR_LDA_OvR.m
                fprintf('Step 3. Individual Electrodes, One vs. Rest\n');
            all_class_sub_indiv_elect_cv_mRMR_LDA_OvR_data.(class_label) = indiv_electrode_cv_mRMR_LDA_OvR(p, class_label);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% Step 4. indiv_electrode_cv_mRMR_LDA.m
                fprintf('Step 4. Individual Electrodes, Regular Labels\n');
            all_class_sub_indiv_elect_cv_mRMR_LDA_data.(class_label) = indiv_electrode_cv_mRMR_LDA(p, class_label);

        end

        %%%% This is a useful block to have, next time you need it it won't take 15
        %%%% years to run. Just load final_checkpoint.mat
        fprintf('Saving checkpoint file... ');
        save(checkpoint_filename_full,...
            'all_class_sub_cv_mRMR_LDA_OvR_data',...
            'all_class_sub_cv_mRMR_LDA_data',...
            'all_class_sub_indiv_elect_cv_mRMR_LDA_OvR_data',...
            'all_class_sub_indiv_elect_cv_mRMR_LDA_data',... 
            '-v7.3');
        fprintf('Done.\n');

    else
        fprintf('Steps 1-4 already saved for %s.  Loading checkpoint file... ', p.current_group_value);
        load(checkpoint_filename_full,...
            'all_class_sub_cv_mRMR_LDA_OvR_data',...
            'all_class_sub_cv_mRMR_LDA_data',...
            'all_class_sub_indiv_elect_cv_mRMR_LDA_OvR_data',...
            'all_class_sub_indiv_elect_cv_mRMR_LDA_data');

        fprintf('Done.\n');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%   Compiling Results   %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Step 5. indiv_electrode_ovr_auc_compiled.m
    fprintf('Step 5. Calculating mean AUC scores and compiling data... ');
    indiv_electrode_ovr_auc_compiled; 
    fprintf('Done.\n');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Step 6. (multinodal scripts)
    %%%%%%%% Only run this if you have to:
    %%%% WARNING: the following function, depending on the number of nodes, may
    %%%%%%%% require multiple hours/days to complete.
    p.auc_results_path = fullfile(p.grouping_path, p.current_group_value, 'AUC_Results');
    auc_files = dir(fullfile(p.auc_results_path, '*.mat'));

    if size(auc_files,  1) ~= p.number_of_nodes
        minutes_elapsed = 0;
        fprintf('Step 6. Now executing multinodal processing to create random data... ');

        job = execute_multinode_pipeline(p);

        while size(auc_files, 1) < p.number_of_nodes + 1
            fprintf('waiting for multinodal analysis...\t%d mins elapsed.\n',  minutes_elapsed);
            minutes_elapsed = minutes_elapsed + 1;
            auc_files = dir(fullfile(p.auc_results_path, '*.mat'));
            if size(auc_files, 1) == p.number_of_nodes
                break
            end
            pause(60)
        end
        fprintf('Done.\nNow compiling results... ');

    else
        fprintf('Step 6. Multinodal analysis already done.\nNow loading and compiling results... ');
    end

    compile_all_nodes_all_iterations;
    fprintf('Done.\n\n');

    %%%% Thresholding real data with respect to significance indicated by
    %%%% permutation test results
    threshold_table = find_significant_auc_threshold(stim_all_nodes_all_iters, onset_all_nodes_all_iters, p.p_value);
    threshd_auc_data = threshold_indiv_data(sub_stim_e_id_auc_all_labels, sub_onset_e_id_auc_all_labels, threshold_table);

    close all; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%   Results and Stats   %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figures_path = fullfile(p.grouping_path, p.current_group_value, 'Figures');
    mkdir(figures_path);

    %%%% Step 7. Visualizing results for all electrodes, regular labels
    all_reg_visuals;        % fig 1
    %%%% Step 8. Visualizing results for all electrodes, ovr labels
    all_ovr_visuals;        % fig 2-4
    %%%% Step 9. Visualizing distribution of aucs based on random data
    perm_test_norm_visual;  % fig 5
    %%%% Step 10. Visualizing step 8, but significant electrodes only
    sig_ovr_visuals;        % fig 6-8
    %%%% Step 11. Visualizing proportions of electrodes considered significant
    ind_ovr_e_proportion;   % fig 9-11
    %%%% Step 12. Projecting significant electrodes onto brain
    % mni_electrode_visual;   % fig 12 (not a figure, creates a surf_show plot)

    close all; 

end

tEnd = seconds(toc(tStart)); 
tEnd.Format = 'hh:mm:ss';
fprintf('Finished.\nTotal elapsed time: %s\n\n', tEnd);

save(sprintf('%s_total_elapsed_time.mat', p.grouping_variable), 'tEnd');