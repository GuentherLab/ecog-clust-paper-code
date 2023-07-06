%% note: this script cannot be used until trialwise single-electrode-exclusion data is found


%%% use bootstrap test to determine whether specific electrodes preferentially encode consonant vs. vowel vs. syllable
%
% Analysis Outline:
% % % for each subject, create n_subsets subsets of trials each containing n_trials_per_subset
% % %   -for each subset: 
% % %       -for each phonemic unit (consonant, vowel, syllable)
% % %       -compute ranks of all electrodes (within this subject), ranging from 0 (weakest) to 1 (strongest), according to the accuracy on this phonunit in this trial subset
% % %   -for each electrode (i.e. when leaving that electrode out):
% % %       -for each subset, compute the difference in ranks for each phonunit pair for this left-out electrode
% % %       -for each subset, compute a "Phonemic-unit Preference Index" (PPI), which is the maximum of (this_phonunit_rank - other_phonunit_rank_1) and (this_phonunit_rank - other_phonunit_rank_2) 
% % %       -(PPI greater than zero indicates that this electrode more strongly encodes this phonunit than both of the other two phonunits)
% % %       -compute 95% CIs for PPI across all subsets; if for any phonunit thesee CIs  do not overlap with zero, then this electrode significantly preferentially encodes this phonunit more than the other 2 phonunits
% % %   -create a brain plot (L hem and R hem) showing electrodes which preferentially encode any of the phonunits (each electrode can preferentially encode at most 1 phonunit), with different colors for each phonunit
%
% note that setting ops.phon_pref_ind_mode='greater_than_one' will allow an electrode to have a nonzero PPI for multiple phonunits
% .... to force electrodes to have nonzero PPI for no more than 1 phonunit, set ops.phon_pref_ind_mode='greater_than_all'
%
% potential usage: if we want to address preferential encoding in labeled areas (ie IFG, STG, M1, etc): at the stage of computing accuracy for a subset, 
%   ... take the mean accuracy for all left-out electrodes in a given labeled area; then instead of assigning ranks for all electrodes, assign ranks to all labeled areas for this subset
%
% because actual trialwise data was not available when writing this ...
% ...... a function gen_dummy_data_for_bootstrap_test.m was created for simulating a simplified version of the dataset
%
% to prepare real trialwise data to be loaded by this script, first run the script: save_trialwise_accuracy_data.m

clear
set_paths();

%%%%%%%%%%% the following option will determine which version of the machine learning analysis data to load
ops.analysis_window_version = '2sec';
% ops.analysis_window_version = 'cluster-specific';

 ops.n_subsets = 20; 
 ops.subset_size_fraction = 1.0; % make each subset this fraction of the size of the total number of trials for this subject.... Alfonso NC recommends using fraction = 1
 ops.replacement_in_subsets = 1; % if true, trials can be used multiple times within a subset; Alfonso NC recommends doing this replacement
 %
 ops.phon_pref_ind_mode = 'greater_than_one'; %% PPI will equal the greatest advantage this phonunit has over all other phonunits (min = zero)
 % ops.phon_pref_ind_mode = 'greater_than_all'; %% PPI will equal the minimum advantage this phonunit has over all other phonunits (min = zero)
 %
 ops.generate_dummy_data = 0;
    genops.all_elcs_acc_boost = 0.1; % add this accuracy when no electrodes are excluded
    genops.elc_table_filename = ALL_ELC_DATA_FILE; 
    genops.acc_offset = -0.3; % increase mean accuracy (from 0.5) by this amount
    genops.ntrials_by_sub = table([357; 362; 369; 372; 376],[108; 143; 118; 144; 133],'VariableNames',{'sub','ntrials'}); % actual ntrials may be lower
    genops.phonunit_names = {'cons','vow','syl'}; 

%
if  ops.generate_dummy_data     % create mock data
    [elc, subdata] = gen_dummy_data_for_bootstrap_test(genops);
elseif ~ops.generate_dummy_data         % load real data
    switch ops.analysis_window_version
        case '2sec'
            load(TRIAL_ACC_FILENAME_2SEC_WINDOW, 'elc', 'subdata')
        case 'cluster-specific'
            load(TRIAL_ACC_FILENAME_CLUSTER_SPEC_WINDOW, 'elc', 'subdata')
    end
end


%%
n_phonunits = length(genops.phonunit_names); 
nsubs = height(subdata);
subdata.n_trials_per_subset = round(subdata.ntrials * ops.subset_size_fraction); 
nelc = height(elc);
for thisphon = genops.phonunit_names
    elc.([thisphon{1},'_pref_sgn']) = false(nelc,1); % binary - whether the 95CI is above zero for this phonunit
    elc.([thisphon{1},'_pref_ind_lower_95CI']) = nan(nelc,1); % PPI lower 95% confidence interval 
    elc.([thisphon{1}, '_subset_acc_change']) = repmat({nan(1,ops.n_subsets)}, nelc, 1);
    elc.([thisphon{1}, '_within_subset_rank']) = nan(nelc,ops.n_subsets); 
    elc.([thisphon{1},'_pref_ind']) = nan(nelc,ops.n_subsets);   % 'Phonemic-unit Preference Index'.... see Analysis Outline
end

%% generate trial subsets, compute accuracy w/out each electrode for each subset
subsets = table(cell(nsubs,ops.n_subsets),'VariableNames',{'trialinds'},'RowNames',cellstr(num2str(subdata.sub)));

for isub = 1:nsubs
    thissub = subdata.sub(isub);
    for isubset = 1:ops.n_subsets
        ntrials_this_sub = subdata.ntrials(isub);
        % shuffled trial indices for this subset
        if  ops.replacement_in_subsets
             trials_this_subset = randi(ntrials_this_sub, 1, subdata.n_trials_per_subset(isub)); 
        elseif ~ops.replacement_in_subsets
            trials_this_subset = randperm(ntrials_this_sub, subdata.n_trials_per_subset(isub)); 
        end
        subsets{num2str(thissub),1}{isubset} = trials_this_subset; % store subset trial indices
        subjrowmatch = find(elc.subject == thissub); 
         nelcs_thissub = length(subjrowmatch); 
         for ielc_within_sub = 1:nelcs_thissub
             irow = subjrowmatch(ielc_within_sub);
             % get accuracy change for this electrode for each phonunit in this subset
             for iphonunit = 1:n_phonunits
                 thisphon = genops.phonunit_names{iphonunit};
                 % find accuracy when leaving this electrode out electrode within this subset of trials for this phonunit
                 subset_acc_without_elc = mean( elc{irow,[thisphon,'_acc_without_elc']}{:}(trials_this_subset) );
                 % find accuracy when all of the subject's electrodes are included within this subset of trials for this phonunit
                 subset_acc_all_elcs = mean( subdata{isub,[thisphon,'_acc_all_elcs']}{:}(trials_this_subset) );
                 % subtract leave-one-out acc from all-elcs acc
                elc{irow,[thisphon,'_subset_acc_change']}{:}(1,isubset) = subset_acc_all_elcs - subset_acc_without_elc; 
             end
         end

         %% convert accuracy change into within-phonunit rankings, so that across-phonunit comparisons are not biased
         %%%% we first make a copy of electrode accuracies for a single subject, then determine within-subject within-subset rankings...
         %%%% ...... then add these values back into the all-subject electrode table
         clear acc_change_this_sub_this_phonunit rankvals_to_assign
         elcs_ordered_within_subset = nan(nelcs_thissub,1); 
         rankvals_to_assign = linspace(0,1,nelcs_thissub)'; 
         for iphonunit = 1:n_phonunits    
            thisphon = genops.phonunit_names{iphonunit};
            acc_change_this_sub_this_phonunit = cell2mat(elc{subjrowmatch, [thisphon,'_subset_acc_change']}); 
            %%%% 'sort' requires 2 steps to get within-subset electrode rankings
            % the worst-ranking electrodes will be listed first within elcs_ordered_within_subset, because we use 'ascend'
            [~, elcs_ordered_within_subset] = sort(acc_change_this_sub_this_phonunit(:,isubset),'ascend');
            % first, get within-subject rankings using rows relative to the list of electrodes within this subject
            ranks_within_subset_subj_relative = nan(nelcs_thissub,1);
            ranks_within_subset_subj_relative(elcs_ordered_within_subset) = rankvals_to_assign; % use elcs_ordered_within_subset to assign 0-1 ranks
            % fill in the all-subject table with subject-relative ranks
             elc{:,[thisphon, '_within_subset_rank']}(subjrowmatch,isubset) = ranks_within_subset_subj_relative; 
         end
    end         

end


 %% compare phonunits within subsets (consonant vs. vowel vs. syllable), get 95% CI for each electrode/phonunit
 % compute 'Phonemic Unit Preference Index' (_pref_ind).... see Analysis Outline
 % ..... compute all electrodes at once
 zeromat = zeros(nelc, ops.n_subsets); %% matrix of zeros to compare against the 2 pairwise advantages
 pairwise_advantage_1 = nan(nelc, ops.n_subsets);
 pairwise_advantage_2 = nan(nelc, ops.n_subsets);
for iphonunit = 1:n_phonunits
     thisphon = genops.phonunit_names{iphonunit};
     other_phon_1 = genops.phonunit_names{mod(iphonunit,n_phonunits)+1};
     other_phon_2 = genops.phonunit_names{mod(iphonunit+1,n_phonunits)+1};
     pairwise_advantage_1 = elc.([thisphon,'_within_subset_rank']) - elc.([other_phon_1,'_within_subset_rank']);
     pairwise_advantage_2 = elc.([thisphon,'_within_subset_rank']) - elc.([other_phon_2,'_within_subset_rank']);
     mat3d = cat(3,zeromat,pairwise_advantage_1,pairwise_advantage_2); % concat to use max()
     
     switch ops.phon_pref_ind_mode
         case 'greater_than_one' % pref_ind will be nonzero if within-subset rank for this electrode is greater for this phonunit than either of the other phonunits
              elc.([thisphon,'_pref_ind']) = max(mat3d,[],3);
         case 'greater_than_all'  % pref_ind will be nonzero if within-subset rank for this electrode is greater for this phonunit than both other phonunits
              elc.([thisphon,'_pref_ind']) = max(zeromat, min(pairwise_advantage_1, pairwise_advantage_2)); 
     end

    %%%%% get get 95% CI for all electrodes for this phonunit
    % the following command sorts PPI values for this phoneme across all subsets (columns) within each electrode (rows)
    pref_ind_w_subsets_sorted =  cell2mat(cellfun(@sort,mat2cell(elc.([thisphon,'_pref_ind']),ones(1,nelc), ops.n_subsets),'UniformOutput',false));
    lower_95CI_subset_index = round(ops.n_subsets * 0.05); 
    elc.([thisphon,'_pref_ind_lower_95CI']) = pref_ind_w_subsets_sorted(:,lower_95CI_subset_index); 
    sgnf_elc_rows = 0 < elc.([thisphon,'_pref_ind_lower_95CI']); % electrodes with PPI significantly greater than zero for this phoneme
    elc{sgnf_elc_rows,[thisphon,'_pref_sgn']} = true; 
end









