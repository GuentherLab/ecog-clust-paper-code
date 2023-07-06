%{

This function will perform CV, mRMR on training set, LDA
training on features selected from training set, and tests the model with
the remaining fold of the CV iteration. 

for each fold, it tabulates the mRMR predictors selected, the predictor
scores, the LDA accuracy for each particular label. 

%}

function alignment_avg_acc = test_LDA_fxn(params, class_label)

topN_feat = params.topN_feat_pool;

all_subs_feat_table = params.pooled_electrode_feature_set;

all_subs_stim_onset = [all_subs_feat_table(:,'sub_num'), all_subs_feat_table.(class_label)];

if mean(all_subs_stim_onset.stim_feat{1}.stim_feat.(1)) == 0
    all_subs_stim_onset.stim_feat = []; 
elseif mean(all_subs_stim_onset.onset_feat{1}.onset_feat.(1)) == 0    
    all_subs_stim_onset.onset_feat = [];
end

sub_nums = all_subs_stim_onset.sub_num; 

for sub_num_idx = 1:length(sub_nums)
    sub_num = sub_nums(sub_num_idx);
    sub_str = strcat('S', string(sub_num));
    
    sub_stim_onset = all_subs_stim_onset(sub_num_idx, :);

    alignment_data = sub_stim_onset.(2){1}; 
    alignment_feat = alignment_data.(1); 
    alignment_class_labels = alignment_data.class_labels; 

    rng(1);
    cvp = cvpartition(alignment_class_labels.(1), 'KFold', 5); 
    
    alignment_label_table = table([1:1:5]', 'VariableNames', {'fold'});
    
    alignment_fold_table = table();

    for fold_num = 1:1:5
        train_idxs = training(cvp, fold_num);
        test_idxs = test(cvp, fold_num);

        X_alignment_train = alignment_feat(train_idxs, :);
        X_alignment_test = alignment_feat(test_idxs, :);

        y_train = alignment_class_labels(train_idxs, :);
        y_test = alignment_class_labels(test_idxs, :);

        alignment_train = [X_alignment_train, y_train];

        [alignment_feat_idx, alignment_feat_scores] = fscmrmr(alignment_train, class_label);

        alignment_feat_names = X_alignment_train.Properties.VariableNames;

        alignment_topN_idx = alignment_feat_idx(1:topN_feat);
        alignment_topN_scores = alignment_feat_scores(alignment_topN_idx);
        alignment_topN_feat_names = alignment_feat_names(alignment_topN_idx);

        alignment_topN_train = X_alignment_train(:, alignment_topN_feat_names);

        alignment_mdl = fitcdiscr(alignment_topN_train, y_train, 'DiscrimType', 'linear', 'ScoreTransform', 'logit');

        [y_pred_alignment, y_pred_probs_alignment, ~] = predict(alignment_mdl, X_alignment_test); 

        alignment_acc = mean(strcmp(y_pred_alignment, y_test.(1)));

        alignment_fold_table_temp = table({alignment_topN_idx}, {alignment_topN_scores}, {alignment_topN_feat_names}, {alignment_mdl}, {y_pred_probs_alignment}, {y_pred_alignment}, {y_test.(1)}, alignment_acc,...
            'VariableNames', {'alignment_topN_idx', 'alignment_topN_scores', 'alignment_topN_feat_names', 'alignment_mdl', 'alignment_y_pred_probs', 'alignment_y_pred', 'alignment_y_test', 'alignment_acc'});
        alignment_fold_table = [alignment_fold_table; alignment_fold_table_temp];

    end
    
    alignment_avg_acc = mean(alignment_fold_table.alignment_acc);
    
end

end
