%%% main loop for creating spatiotemporal matched filters and using them to classify the feature val for each trial
% 
% called by stim_classification.m
%
% updated 2021/5/6

% create STMFs for each feature val
for iftrval = 1:n_feature_vals %%% iftrval = index of whatever feature is being analyzed - full word, consonant pair, vowel, etc. 
    match_trials = trials{:, analysis_ops.feature_to_analyze} == iftrval; %%% trials for this subject across all blocks which have this feature value
    ecog_match_trials = trials.ecog_trial_ind(match_trials); % get the matching index within preprocessed_data
    %%%% stmf = mean response to all trials with this feature val, for each electrode
    %%% use sum rather than mean because it will make subtracting the left-out trial simpler
    feature_vals.stmf{iftrval} = sum( preprocessed_data.data(sites_to_analyze,...
        analysis_ops.window_start_stop(1):analysis_ops.window_start_stop(2), ecog_match_trials), 3); 
end

ftrval_r_cor = NaN(ntrials, n_feature_vals); % pearson correlation coefficient between each trial and STMFs
best_iftrval = NaN(ntrials, 1); 
for itrial = 1:ntrials
    correct_ftrval_ind = trials{itrial, analysis_ops.feature_to_analyze}; % get the feature val index present in this trial
    for iftrval = 1:n_feature_vals
        % subtract the current trial's responses from the STMF which it was used to construct (leave-one-out classification)
        %%% to save time, don't divide by ntrials to get mean because it won't change the correlation
        ecog_trial_ind = trials.ecog_trial_ind(itrial); % index of this trial within preprocessed_data
        if  iftrval == correct_ftrval_ind % if this trial was used to construct this STMF
            stmf_leave_one_out = feature_vals.stmf{iftrval} - ...
                preprocessed_data.data(sites_to_analyze, analysis_ops.window_start_stop(1):analysis_ops.window_start_stop(2), ecog_trial_ind); 
        elseif iftrval ~= correct_ftrval_ind%% if this trial wasn't included in this STMF
            stmf_leave_one_out = feature_vals.stmf{iftrval} ; 
        end
        cor = corrcoef( stmf_leave_one_out, ...  % get cor
            preprocessed_data.data(sites_to_analyze,analysis_ops.window_start_stop(1):analysis_ops.window_start_stop(2),ecog_trial_ind) );
        ftrval_r_cor(itrial, iftrval) = cor(2,1);
        [~, best_iftrval(itrial)] = max(ftrval_r_cor(itrial,:));
    end
end
trials.correct_ftrval(:, iter) = trials{:, analysis_ops.feature_to_analyze} == best_iftrval;  % check whether classified ftrval matches actual ftrval