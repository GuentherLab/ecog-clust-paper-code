%%%% classify stim features using spatio-temporal matched filters
%
% call this function from stmf_toplevel with ops structure (see stmf_toplevel for ops input)
%%% results will be saved into savedir as .mat files, not output into workspace
% save_append will be added the output filename
%
% methods adapted from Ramsey et al. 2018 (pubmed 28993231)
% 
%%% updated by AM 2021/5/27

function stmf_classification(analysis_ops, savedir, save_append)

%% set defaults
vardefault('savedir', '/projectnb/busplab/Experiments/ECoG_Preprocessed_AM/stmf_results'); % path to save results in
vardefault('save_append', datestr(clock)); 
vardefault('analysis_ops',struct); 

field_default('analysis_ops',  'n_iterations', 1)  %% run this many iterations on each node
field_default('analysis_ops', 'subjects_to_analyze', [357, 362, 369, 372, 376]); 

field_default('analysis_ops',  'feature_to_analyze', 'word')

field_default('analysis_ops',  'alignment_condition', 'onset')
% field_default('analysis_ops',  'alignment_condition', 'stimuli')

%%%%% window edges in samples to analyze from each trial
%%%% Ramsey ea 2018 used -500 to 500ms peri-voice-onset
field_default('analysis_ops', 'window_start_stop', [1500 2500])

 % analyze_full_electrode_set==1: in addition to analyzing subsets of electrodes, do an analysis using all electrodes
 %%% this full electrode set will NOT include electrodes excluded by only_analyze_clustered_sites and use_both_align_conditions
field_default('analysis_ops', 'analyze_full_electrode_set', 0)
field_default('analysis_ops','only_analyze_clustered_sites', 1) % for classification, only use sites that were assigned a cluster
    field_default('analysis_ops', 'use_both_align_conditions', 0) % has effect only if only_analyze_clustered_sites == true
    
%%% paths
% folder with analyzed ecog data for each sub
field_default('analysis_ops', 'ecog_localepoched_topdir', '/projectnb/busplab/Experiments/ECoG_Preprocessed/LocalEpoched')
     % .mat file containing ecog responses
    field_default('analysis_ops', 'ecog_data_filename', ['Epoch_', analysis_ops.alignment_condition, '_1_Hilbert_HG'])
%%% path to AnalysisClass.m, which we need to in order to load ecog data
field_default('analysis_ops', 'analysis_class_path', '/project/busplab/software/ecog')
field_default('analysis_ops', 'electrode_list_file', '/projectnb/busplab/Experiments/ECoG_Preprocessed_AM/electrodes.mat')
field_default('analysis_ops', 'clusterkey_file', '/projectnb/busplab/Experiments/ECoG_Preprocessed_AM/clusterkey')

    

%% analysis
starttime = datestr(clock); 
addpath(analysis_ops.analysis_class_path)
match_trial_2_behavior()
nsubs = length(analysis_ops.subjects_to_analyze); 
stmf_results = table(analysis_ops.subjects_to_analyze',cell(nsubs,1),cell(nsubs,1), 'VariableNames',... % results for all subjects
                                                               {'subject',       'trials',           'sites'}); 
rng('shuffle'); % make sure we don't repeat seeds when running on a new matlab session

for isub = 1:nsubs
    this_sub = analysis_ops.subjects_to_analyze(isub);
    clear preprocessed_data feature_vals trials iter_results
    load([analysis_ops.ecog_localepoched_topdir, filesep, 'S', num2str(this_sub), filesep, analysis_ops.ecog_data_filename])   %%% load ecog data
    sample_rate_hz = preprocessed_data.Fs; 
    feature_vals = stim_index_key(strcmp(stim_index_key.type, analysis_ops.feature_to_analyze), :); %%% all values which this stim feature can take on
    n_feature_vals = height(feature_vals); 
    feature_vals.stmf = cell(n_feature_vals,1); 
    trials = filelist.trials{filelist.subject==this_sub}; % stim data for each trial we're about to analyze from this subject
    iter_results = table(NaN(analysis_ops.n_iterations,1), cell(analysis_ops.n_iterations,1), NaN(analysis_ops.n_iterations,1), ...
        'VariableNames', {'nsites', 'idx_LocalProcessed', 'frac_correct'}); 

    % check that the number of trials from our excel file list matches the number of ecog trials in Local Epoched folder
    if height(trials) ~= size(preprocessed_data.data, 3)  % height(trials) is the number of trials listed in excel file list
        error('number of trials in excel file list does not match number of trials in Local Processed ecog data folder')
    end
    
    trials.ecog_trial_ind = [1:height(trials)]'; % save the index for matching to ecog data before discarding bad trials
    trials = trials(trials.analyze, :); % discard trials which were marked as bad
    ntrials = height(trials);
    
    %% determine which sites to use
    %%% values in 'available_sites' are rows of preprocessed_data, NOT the chan labels in preprocessed_data.chan_ids
    %%% available_sites lists all sites which may be used for classification...
    %%% ...a smaller subset (sites_to_analyze) may be used on individual iterations
    if strcmp(analysis_ops.alignment_condition, 'onset')
        align_idx = 1; % matches variable 'type' in electrodes table
    elseif strcmp(analysis_ops.alignment_condition, 'stimuli')
        align_idx = 2; % matches variable 'type' in electrodes table
    end
    
    load(analysis_ops.electrode_list_file, 'electrodes')
    matchrows = electrodes.subject == this_sub; % find sites in electrodes table which match this subject
    if ~analysis_ops.only_analyze_clustered_sites
        available_sites = 1:preprocessed_data.num_chans; % use all sites, including unclustered
    elseif analysis_ops.only_analyze_clustered_sites
        matchrows = matchrows  & ~isnan(electrodes.cluster);  % find clustered sites which match this subject
        
         if ~analysis_ops.use_both_align_conditions % if we're only using 1 align condition for classification
            matchrows = matchrows & electrodes.type == align_idx; % only include this align condition
         end
        
        available_sites = electrodes.idx_LocalProcessed(matchrows); % clustered sites.... see match_clustered_to_LocalEpoched.m
        available_sites = unique(available_sites); % sort, remove double-listed sites (for 2 conditions)
    end
    nsites_max = length(available_sites); 
    
    % prepare table of sites with cluster labels
    sites = table; 
    for isite = 1:nsites_max
        thissite = available_sites(isite); % index in LocalProcessed
        sitematch = matchrows & electrodes.idx_LocalProcessed == thissite; % find row in electrodes table
        % if this site was listed under both conditions and we're analyzing sites clustered under both conditions...
        % .... use only the site entry corresponding to our favored align condition
        if nnz(sitematch) > 1 % if there's a double entry
            sitematch = sitematch & electrodes.type == align_idx; % match align condition
        end
        
        % if this site is not listed in the electrodes table, give it an empty row
        %%% if it is listed, add its row from electrodes to our new sites table
        if nnz(sitematch) == 0
            emptyrow = table(this_sub,        thissite,                  {''},                false,                      NaN,        NaN, align_idx, {''},       {''},           2,             {''},...
                'VariableNames', ...
                                    {'subject', 'idx_LocalProcessed', 'cluster_name', 'multiple_clusters', 'electrode', 'cluster', 'type',    'roi', 'region', 'right_hemi', 'notes'});
            sites = [sites; emptyrow]; % append
        elseif nnz(sitematch)  > 0 % this site is listed in electrodes table
            sites = [sites; electrodes(sitematch,:)]; % append
        end
        
    end
    
    % make clearer variable labeling hemisphere
    sites.hemi = cell(nsites_max, 1);
    sites.hemi(sites.right_hemi == 1) = {'r'};
    sites.hemi(sites.right_hemi == 2) = {''}; % this site was not clustered, therefore hemi not yet labeled
    sites.hemi(isnan(sites.right_hemi)) = {'l'};
    sites.right_hemi = []; % delete old var
    
    % cluster labeling; we will use the clusterkey.unv_clust (universal cluster) labels, not the align-condition-speicfic labels
    load(analysis_ops.clusterkey_file, 'clusterkey')
    sites.cluster = []; % delete the old condition-specific version of cluster indices
    sites.clust = NaN(nsites_max, 1); % universal cluster labels
    for iclust = 1:height(clusterkey)
        namematch =  strcmp( sites.cluster_name, clusterkey.name(iclust));
        sites.clust(namematch) = clusterkey.unv_clust(iclust); 
    end
    sites = movevars(sites, {'clust', 'cluster_name'}, 'Before', 'idx_LocalProcessed');
    stmf_results.sites{isub} = sites; 
    
    %% run stmf classification with all sites
    if analysis_ops.analyze_full_electrode_set
        sites_to_analyze = available_sites;   
        trials_orig = trials; 
        trials.correct_ftrval = false(ntrials, 1); % whether the STMFs correctly classified the feature val in this trial
        iter = 1; % use just 1 iteration for full electrode set
        stmf_loop(); % run once
        trials_allsites = trials; % rename
        trials = trials_orig; 
        trials_allsites.ftrval_r_cor = ftrval_r_cor; % save results in table
        trials_allsites.best_iftrval = best_iftrval;  
        stmf_results.trials_allsites{isub} = trials_allsites; % save into results table
        clear trials_orig sites_to_analyze ftrval_r_cor best_iftrval trials_allsites
    end
    
    %% run classification for all electrode subsets
    trials.correct_ftrval = false(ntrials, analysis_ops.n_iterations); % whether the STMFs correcly classified the feature val in this trial
    for iter = 1:analysis_ops.n_iterations
        n_sites_this_iter = randperm(nsites_max,1); % randomize the number of electrodes to use for this iteration
        avail_inds = randperm(nsites_max, n_sites_this_iter); % random site indices (relative to available_sites)
        sites_to_analyze = available_sites(avail_inds); % get electrodes subset for this iter; site indices from preprocessed_data.data
        stmf_loop(); 
        iter_results.nsites(iter) = n_sites_this_iter;
        iter_results.idx_LocalProcessed{iter} = sites_to_analyze; 
        iter_results.frac_correct(iter) = mean(trials.correct_ftrval(:, iter)); % get fraction correct classification from this iteration
    end
    stmf_results.trials{isub} = trials; % add results to table
    stmf_results.iter_results{isub} = iter_results; 
    % mean(trials.correct_ftrval) % display results
end

% save results
vardefault('save_append', ''); 
field_default('analysis_ops', 'save_prepend', 'stmf_results_'); 
savename = [savedir, filesep, analysis_ops.save_prepend, analysis_ops.feature_to_analyze, save_append]; 
endtime = datestr(clock); 
save(savename, 'stmf_results', 'analysis_ops', 'starttime', 'endtime') 
