%%% match channel indices in LocalEpoched data with indices in clustered data...
%%% ... then add MNI coordinates from LocalProcessed to each electrode
%
% electrodes in LocalEpoched were selected based on simple responsiveness...
% ... many of these electrodes get excluded during Kalman filtering/clustering
%
% this script takes the clustered electrodes and matches them with their location in LocalEpoched
% (variable name = preprocessed_data.data) data matrix, which also contains the electrodes that
% were excluded before final cluster labels were assigned
%
  % this script also adds a variable 'global_clust_num' to the electrodes.mat table which corresponds to a given cluster... 
  % ... across both align conditions; numbers will be in order of onset-aligned response start time, ...
  % ... from fig. 2 / table 2 of Scott Kuzdeba's dissertation
  % these cluster numbers can then be called by plot_greenlee_electrodes_on_brain.m to plot cluster label consistently....
  % .... given the same colormap is used across plots
%
% updated by Andrew Meier 2021/8/11

%% set directories

subject_list = [357, 362, 369, 372, 376]; % subjects to get data from

% this electrode list will be appended with site indices from LocalEpoched
% % % % % % % % % % % % % % % % electrodes_excel = '/projectnb/busplab/Experiments/ECoG_Preprocessed_AM/old/electrodes.xlsx'; 
electrodes_table_file = '/projectnb2/busplab/Experiments/ECoG_Preprocessed_AM/electrodes.mat'; 

% chan ids from LocalEpoched (not clustered) will be loaded from here
% ...this file was originally created by Generate_Chan_IDs_mat.m
chan_ids_file =  'project/busplab/software/ecog/data/LocalEpoched_chan_ids'; 

% epoched data (before Kalman filtering and clustering) is stored here
LocalEpoched_dir = '/projectnb/busplab/Experiments/ECoG_Preprocessed/LocalEpoched'; 

% list of cluster labels
clusterkey_file = '/project/busplab/software/ecog/data/clusterkey'; 

% directory containing MNI coordinates for each electrode
LocalProcessed_dir = '/projectnb2/busplab/Experiments/ECoG_Preprocessed/LocalProcessed'; 

%% get cluster labels
% % % % % % % % % % % % % % % % % % electrodes = readtable(electrodes_excel); % read site list with cluster labels
load(electrodes_table_file)
electrodes.idx_LocalProcessed = NaN(height(electrodes), 1); % electrode index within preprocessed_data.data
electrodes.cluster_name = cell(height(electrodes), 1); % cluster name... should be consistent across alignment conditions
electrodes.multiple_clusters = false(height(electrodes), 1); % true if this electrode has different cluster across alignment conditions
electrodes.idx_LocalProcessed = NaN(height(electrodes), 1); % electrode index within preprocessed_data.data
electrodes.global_clust_num = NaN(height(electrodes), 1); % cluster numbers in order of onset-aligned response start time
load(chan_ids_file); % struct will be loaded as 'chan_list' variable
load(clusterkey_file, 'clusterkey')

% match channel indices
nsubs = length(subject_list); 
for isub = 1:nsubs
    this_sub = subject_list(isub); 
    this_sub_name = num2str(this_sub); 
    these_LocalEpoched_chans = chan_list.(['S', this_sub_name]); % chan indices in preprocessed_data.data
    for ichan_Local = 1:length(these_LocalEpoched_chans)
        channame = these_LocalEpoched_chans(ichan_Local); % label of this chan from preprocessed_data.chan_names
        matchrows = electrodes.subject == this_sub & electrodes.electrode == channame; % matching rows in electrodes table
        electrodes.idx_LocalProcessed(matchrows) = ichan_Local; % corresponding index from preprocessed_data.data
    end
end

% add cluster labels
for iclust = 1:height(clusterkey)
    matchrows = electrodes.type == clusterkey.align_type(iclust) & electrodes.cluster == clusterkey.clust_ind(iclust); 
    clustnames = electrodes.cluster_name; % deal was causing problems with table var, use temp var instead
    [clustnames{matchrows}] = deal(clusterkey.name{iclust}); 
    electrodes.cluster_name = clustnames; 
    % add global cluster indices
    globalnums = electrodes.global_clust_num; 
    [globalnums(matchrows)] = deal(clusterkey.global_clust_num(iclust)); 
    electrodes.global_clust_num = globalnums; 
end

% check that all entries for a single electrode have the same cluster label
for irow = 1:height(electrodes)
    matchrows = electrodes.subject == electrodes.subject(irow) & electrodes.electrode == electrodes.electrode(irow); 
    clustnames = electrodes.cluster_name(matchrows); 
    if ~any(cellfun(@isempty, clustnames)) && length(unique(clustnames)) > 1 % if this site has different cluster across conditions
        electrodes.multiple_clusters(matchrows) = true;
    end
end

electrodes = movevars(electrodes, {'idx_LocalProcessed', 'cluster_name', 'multiple_clusters'},...
    'After', 'subject'); % rearrange table

%% get MNI coordinates
electrodes.mni_coord = NaN(height(electrodes),3);
for isub = 1:nsubs
    this_sub = subject_list(isub); 
    this_sub_name = num2str(this_sub); 
    load([LocalProcessed_dir, filesep, 'S', this_sub_name, filesep, 'S', this_sub_name, '_MNI'])
    thissub_rows = find(electrodes.subject == this_sub);
    for grandtable_row = thissub_rows' % iterate over every index in electrodes table with electrode from this subject
        electrode_old_ind = electrodes.electrode(grandtable_row); % index to match with MNI table indexing scheme i.e. MNI(:,1)
        mni_matchrow = electrode_old_ind == MNI(:,1); 
        electrodes.mni_coord(grandtable_row,:) = MNI(mni_matchrow, 2:4); % fill in mni coordinates
    end
end
electrodes = movevars(electrodes, {'mni_coord'}, 'Before', 'roi'); % rearrange table

note_source = 'original electrode data from Scott Kuzdeba stored in ECoG-Clust Google Drive in ''old'' folder under: All_Electrode_labels.xlsx and Electrodes_README.txt';
note_clusters_mni = 'electrode table modifedy with function /project/busplab/software/ecog/classification_AM/match_clustered_and_MNI_to_LocalEpoched.m';
clearvars('-except', 'electrodes', 'note_source', 'note_clusters_mni')
