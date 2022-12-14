%{

This function will perform CV, mRMR on training set, LDA
training on features selected from training set, and tests the model with
the remaining fold of the CV iteration. 

for each fold, it tabulates the mRMR predictors selected, the predictor
scores, the LDA accuracy for each particular label. 

%}

function sub_cv_mRMR_LDA_data = all_electrodes_cv_mRMR_LDA(params, class_label)

cv_mRMR_LDA_res_file = 'sub_cv_mRMR_LDA_data.mat';
cv_mRMR_LDA_res_filename_full = fullfile(params.grouping_path, params.current_group_value, class_label, cv_mRMR_LDA_res_file);

if ~isfile(cv_mRMR_LDA_res_filename_full)
topN_feat = params.topN_feat_pool;

all_subs_feat_table = params.pooled_electrode_feature_set;

all_subs_stim_onset = [all_subs_feat_table(:,'sub_num'), all_subs_feat_table.(class_label)];

sub_nums = all_subs_stim_onset.sub_num; 

sub_cv_mRMR_LDA_data = struct();

for sub_num_idx = 1:length(sub_nums)
    sub_num = sub_nums(sub_num_idx);
    sub_str = strcat('S', string(sub_num));
    fprintf('\tSubject %s\n', sub_str);
    
    sub_stim_onset = all_subs_stim_onset(sub_num_idx, :);
    stim_feat = sub_stim_onset.stim_feat{1}.stim_feat;
    onset_feat = sub_stim_onset.onset_feat{1}.onset_feat;
    
    class_labels = sub_stim_onset.stim_feat{1}.class_labels;

    rng(1);
    cvp = cvpartition(class_labels.(1), 'KFold', 5); 
    
    stim_label_table = table([1:1:5]', 'VariableNames', {'fold'});
    onset_label_table = table([1:1:5]', 'VariableNames', {'fold'});
    
    stim_fold_table = table();
    onset_fold_table = table();

    for fold_num = 1:1:5
        train_idxs = training(cvp, fold_num);
        test_idxs = test(cvp, fold_num);

        X_stim_train = stim_feat(train_idxs, :);
        X_stim_test = stim_feat(test_idxs, :);

        X_onset_train = onset_feat(train_idxs, :);
        X_onset_test = onset_feat(test_idxs, :);

        y_train = class_labels(train_idxs, :);
        y_test = class_labels(test_idxs, :);

        stim_train = [X_stim_train, y_train];
        onset_train = [X_onset_train, y_train];

        [stim_feat_idx, stim_feat_scores] = fscmrmr(stim_train, class_label);
        [onset_feat_idx, onset_feat_scores] = fscmrmr(onset_train, class_label);

        stim_feat_names = X_stim_train.Properties.VariableNames;
        onset_feat_names = X_onset_train.Properties.VariableNames;

        stim_topN_idx = stim_feat_idx(1:topN_feat);
        stim_topN_scores = stim_feat_scores(stim_topN_idx);
        stim_topN_feat_names = stim_feat_names(stim_topN_idx);

        onset_topN_idx = onset_feat_idx(1:topN_feat);
        onset_topN_scores = onset_feat_scores(onset_topN_idx);
        onset_topN_feat_names = onset_feat_names(onset_topN_idx);

        stim_topN_train = X_stim_train(:, stim_topN_feat_names);
        onset_topN_train = X_onset_train(:, onset_topN_feat_names);

        stim_mdl = fitcdiscr(stim_topN_train, y_train, 'DiscrimType', 'linear', 'ScoreTransform', 'logit');
        onset_mdl = fitcdiscr(onset_topN_train, y_train, 'DiscrimType', 'linear', 'ScoreTransform', 'logit');

        [y_pred_stim, y_pred_probs_stim, ~] = predict(stim_mdl, X_stim_test); 
        [y_pred_onset, y_pred_probs_onset, ~] = predict(onset_mdl, X_onset_test);

        stim_acc = mean(strcmp(y_pred_stim, y_test.(1)));
        onset_acc = mean(strcmp(y_pred_onset, y_test.(1)));

        stim_fold_table_temp = table({stim_topN_idx}, {stim_topN_scores}, {stim_topN_feat_names}, {stim_mdl}, {y_pred_probs_stim}, {y_pred_stim}, {y_test.(1)}, stim_acc,...
            'VariableNames', {'stim_topN_idx', 'stim_topN_scores', 'stim_topN_feat_names', 'stim_mdl', 'stim_y_pred_probs', 'stim_y_pred', 'stim_y_test', 'stim_acc'});
        stim_fold_table = [stim_fold_table; stim_fold_table_temp];

        onset_fold_table_temp = table({onset_topN_idx}, {onset_topN_scores}, {onset_topN_feat_names}, {onset_mdl}, {y_pred_probs_onset}, {y_pred_onset}, {y_test.(1)}, onset_acc,...
            'VariableNames', {'onset_topN_idx', 'onset_topN_scores', 'onset_topN_feat_names', 'onset_mdl', 'onset_y_pred_probs', 'onset_y_pred', 'onset_y_test', 'onset_acc'});
        onset_fold_table = [onset_fold_table; onset_fold_table_temp];

    end
    stim_all_folds_class = table(stim_fold_table, 'VariableNames', {class_label});
    stim_label_table = [stim_label_table, stim_all_folds_class];

    onset_all_folds_class = table(onset_fold_table, 'VariableNames', {class_label});
    onset_label_table = [onset_label_table, onset_all_folds_class];

    sub_cv_mRMR_LDA_data.(sub_str).stim_data = stim_label_table;
    sub_cv_mRMR_LDA_data.(sub_str).onset_data = onset_label_table;
    
    fprintf('\t\tStimulus data complete. | ');    
    fprintf('Onset data complete.\n');
end
fprintf('\n\tComplete. Saving... ');
mkdir(fullfile(params.grouping_path, params.current_group_value, class_label));
save(cv_mRMR_LDA_res_filename_full, 'sub_cv_mRMR_LDA_data', '-v7.3');    
fprintf('Done.\n');    

else
fprintf('\n\Loading sub_cv_mRMR_LDA_data...');
    
sub_cv_mRMR_LDA_data = load(cv_mRMR_LDA_res_filename_full).sub_cv_mRMR_LDA_data;    
fprintf('Done.\n');
end

end
