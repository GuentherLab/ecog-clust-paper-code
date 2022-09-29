%%% run stmf_classification for many iterations on a computing cluster
%
% this function calls stmf_classification.m to run many iterations on multiple computing cluster nodes
%%% total number of stmf iterations will be n_nodes*ops.n_iterations
%
% use conn toolbox to submit job to cluster
%
% updated by Andrew Meier 2021/5/16

%% parameters
homedir = '/usr2/postdoc/amsmeier'; %% conn sometimes generates errors if not called from starting directory

n_nodes = 50; %%% run the function on this many computing cluster nodes
% features_to_analyze = {'word', 'consonants', 'vowel'}; 
    features_to_analyze = {'vowel'}; 
ops.n_iterations = 100; %% run this many iterations on each node
ops.subjects_to_analyze = [357, 362, 369, 372, 376]; 

ops.alignment_condition = 'onset'; %%% response window -2 to 1sec around speech onset
%     ops.alignment_condition = 'stimuli'; %%% response window -1 to 2sec around GO cue

%%%%% window edges in samples to analyze from each trial
%%%% Ramsey ea 2018 used -500 to 500ms peri-voice-onset
% ops.window_start_stop = [1500 2500]; % for 'onset' condition, these values replicate the main window from Ramsey ea 2018
ops.window_start_stop = [1000 3000]; 
% ops.window_start_stop = [1 3000]; % use full 3sec window of data

savedir = '/projectnb/busplab/Experiments/ECoG_Preprocessed_AM/stmf_results'; % path to save results in
ops.save_prepend = 'stmf_results2sec_'; 

ops.only_analyze_clustered_sites = 1; % for classification, only use sites that were assigned a cluster
    % if only_analyze_clustered_sites==true, this option toggles using for classificatoin sites that were assigned a
    % cluster under either align condition; this option has no effect if nly_analyze_clustered_sites==false
    ops.use_both_align_conditions = 0; 
 % in addition to analyzing subsets of electrodes, do an analysis using all electrodes
 %%% this full electrode set will NOT include electrodes excluded by only_analyze_clustered_sites and use_both_align_conditions
ops.analyze_full_electrode_set = 0; 

% paths
ops.ecog_localepoched_topdir = '/projectnb/busplab/Experiments/ECoG_Preprocessed/LocalEpoched';   % folder with analyzed ecog data for each sub
    ops.ecog_data_filename = ['Epoch_', ops.alignment_condition, '_1_Hilbert_HG']; % .mat file containing ecog responses
ops.analysis_class_path = '/project/busplab/software/ecog'; %%% path to AnalysisClass.m, which we need to in order to load ecog data
ops.electrode_list_file = '/projectnb/busplab/Experiments/ECoG_Preprocessed_AM/electrodes.mat'; 
ops.clusterkey_file = '/projectnb/busplab/Experiments/ECoG_Preprocessed_AM/clusterkey'; 

%% run analysis
start_time = clock; datestr(start_time)
origdir = pwd; 
cd(homedir)
for ifeature = 1:length(features_to_analyze)
    ops.feature_to_analyze = features_to_analyze{ifeature}; 
    for inode = 1:n_nodes
        qlog(inode) = conn('submit', @stmf_classification, ops, savedir, ['_', num2str(inode)]);
    end
    
    [~, qout] = system('qstat -u amsmeier'); % get job statuses; should be empty if no jobs running
    while ~isempty(qout) % if there's still a job running for this feature
        pause(5) % wait, then check again if jobs are done
        [~, qout] = system('qstat -u amsmeier'); % get job statuses; should be empty if no jobs running
    end
    
    feature_to_analyze = ops.feature_to_analyze 
    stmf_compile_results(); % gather results from all output files
end
cd(origdir); clear origdir
% end_time = clock
