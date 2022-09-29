%%% run this script to organize data in compact form for surf plotting clusters
%
% run this script before running surf_top_coders_per_cluster and plot_right_vs_left
%
% updated by AM 2022/7/5


set_paths()

top_proportion_electrodes = 0.3; % use this proportion of the electrodes as 'top coders'
preloaded_electrode_data_filename = [ROOT_DIR, filesep, 'projectnb/busplab/Experiments/ECoG_Preprocessed_AM/leave_one_out_data'];
load([ROOT_DIR, filesep, 'project/busplab/software/ecog/data/clusterkey.mat'], 'clusterkey')
    global_clust_num_list = [3 4 5 6 7 8]; % voice onset-aligned clusters only
    
    load(preloaded_electrode_data_filename,'ow_sorted') % onset data sorted by word coding
elc = sortrows(ow_sorted,{'subject','electrode'}); % sort electrodes by subject and electrode ID
elc.cluster_name = string(elc.cluster_name);

savename = [ROOT_DIR, filesep, 'projectnb/busplab/Experiments/ECoG_Preprocessed_AM/topelc_data_to_surf']; 
    vars_to_save = {'elc', 'topelc', 'top_proportion_electrodes', 'top_proportion_electrodes', 'clusterkey','global_clust_num_list'};

%%
 %%%%%%%%%%%%%%%% organize data %%%%%%%%%%
n_elcs = height(elc);  
falsecol = false(n_elcs,1);
elc = [elc, table(falsecol, falsecol, falsecol, falsecol, 'VariableNames', {'top_cons_coder','top_vow_coder','top_word_coder','top_cons_or_word'})];                              
elc.feat_set = [];

% create table with top-30% coding electrodes data by cluster
clustlist = {'PtM-s','PtM-r','ME-sb','ME-sn','AP-r','AP-s'}';
nclusts = length(clustlist); 

% make tables of top coders
top_inds = round(n_elcs*[1-top_proportion_electrodes]) : n_elcs; 
n_top_elc = length(top_inds);
elc = sortrows(elc,'cons_accuracy_change_wo_electrode');
    elc.top_cons_coder(top_inds) = true; 
    elc_cons = elc(top_inds,:); 
elc = sortrows(elc,'vowel_accuracy_change_wo_electrode');
    elc.top_vow_coder(top_inds) = true; 
    elc_vow = elc(top_inds,:); 
elc = sortrows(elc,'word_accuracy_change_wo_electrode');
    elc.top_word_coder(top_inds) = true; 
    elc_word = elc(top_inds,:); 
    
elc.top_cons_or_word(elc.top_cons_coder | elc.top_word_coder) = true; 
    elc_cons_or_word = elc(elc.top_cons_or_word,:); 
    elc = sortrows(elc,{'subject','electrode'}); % sort electrodes by subject and electrode ID
    
nancol = NaN(nclusts, 1);
nan2 = [nancol, nancol]; 
topelc = table(clustlist, nancol,   nancol,              nan2,       nancol,              nan2,         nancol,            nan2,       nancol,          nancol,  'VariableNames',...
                    {'clust',       'n_elc', 'cons_left_prop', 'cons_RL', 'vow_left_prop', 'vow_RL', 'word_left_prop', 'word_RL', 'cons_word_left_n', 'cons_word_left_prop'}); 
    
for iclust = 1:nclusts
    thisclust = clustlist{iclust};
    clustrows = strcmp(elc.cluster_name, thisclust);
    topelc.n_elc(iclust) = nnz(clustrows);
    
    clustrows_top_cons = clustrows & elc.top_cons_coder; 
    n_top_cons_this_clust = nnz(clustrows_top_cons); % number of electrodes in this cluster which are top cons coderclc
    
    n_top_cons_this_clust_leftside = nnz(clustrows_top_cons & isnan(elc.right_hemi)); % number electrodes in this clust which are top coders and in L hem 
    topelc.cons_RL(iclust,:) = [n_top_cons_this_clust-n_top_cons_this_clust_leftside  , n_top_cons_this_clust_leftside]; % number of R and L electrodes
    topelc.cons_left_prop(iclust) = n_top_cons_this_clust_leftside / n_top_cons_this_clust; 
    
    clustrows_top_vow = clustrows & elc.top_vow_coder; 
    n_top_vow_this_clust = nnz(clustrows_top_vow); % number of electrodes in this cluster which are top coder
    n_top_vow_this_clust_leftside = nnz(clustrows_top_vow & isnan(elc.right_hemi)); % number electrodes in this clust which are top coders and in L hem 
    topelc.vow_RL(iclust,:) = [n_top_vow_this_clust-n_top_vow_this_clust_leftside  , n_top_vow_this_clust_leftside]; % number of R and L electrodes
    topelc.vow_left_prop(iclust) = n_top_vow_this_clust_leftside / n_top_vow_this_clust; 
    
    clustrows_top_word = clustrows & elc.top_word_coder; 
    n_top_word_this_clust = nnz(clustrows_top_word); % number of electrodes in this cluster which are top coder
    n_top_word_this_clust_leftside = nnz(clustrows_top_word & isnan(elc.right_hemi));  
    topelc.word_RL(iclust,:) = [n_top_word_this_clust-n_top_word_this_clust_leftside  , n_top_word_this_clust_leftside]; % number of R and L electrodes
    topelc.word_left_prop(iclust) = n_top_word_this_clust_leftside / n_top_word_this_clust; 
    
    clustrows_top_cons_or_word = clustrows & elc.top_cons_or_word; 
    n_top_cons_or_word_this_clust = nnz(clustrows_top_cons_or_word); % number of electrodes in this cluster which are top coder
    topelc.cons_word_n(iclust,1) = n_top_cons_or_word_this_clust; 
    n_top_cons_or_word_this_clust_leftside = nnz(clustrows_top_cons_or_word & isnan(elc.right_hemi));
    topelc.cons_word_left_n(iclust) = n_top_cons_or_word_this_clust_leftside; 
    topelc.cons_word_left_prop(iclust) = n_top_cons_or_word_this_clust_leftside / n_top_cons_or_word_this_clust; 
end

%%
save(savename, vars_to_save{:})
