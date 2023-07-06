% load and organize top-encoding electrodes for consonant, vowel, syllable
% only loads onset-aligned, not stim-aligned data
%
% % % run this script before running surf_top_coders_per_cluster and plot_right_vs_left
%
% updated by AM 2022/8/9

function [elc, topelc, clustlist, n_elcs, nclusts] = load_top_encoders(ops)

set_paths()
vardefault('ops',struct);

field_default('ops','top_proportion_electrodes', 0.3); % use this proportion of the electrodes as 'top coders'
% field_default('ops','top_proportion_electrodes', 0.6); % use this proportion of the electrodes as 'top coders'

field_default('ops','preloaded_electrode_data_filename' , [ROOT_DIR, filesep, 'projectnb/busplab/Experiments/ECoG_Preprocessed_AM/leave_one_out_data']);

load(ops.preloaded_electrode_data_filename,'ow_sorted') % onset data sorted by word coding
elc = sortrows(ow_sorted,{'subject','electrode'}); % sort electrodes by subject and electrode ID
elc.cluster_name = string(elc.cluster_name);
load([ROOT_DIR, filesep, 'project/busplab/software/ecog/ecog_clust/data/clust_timing_kuzdeba_table_2.mat'], 'clusterkey')
    global_clust_num_list = [3 4 5 6 7 8]; % voice onset-aligned clusters only
    
savename = [ROOT_DIR, filesep, 'projectnb/busplab/Experiments/ECoG_Preprocessed_AM/topelc_data_to_surf']; 
    vars_to_save = {'elc', 'elc_cons', 'elc_vow', 'elc_word', 'elc_anytop', 'topelc', 'top_proportion_electrodes',...
        'clustlist', 'clusterkey','global_clust_num_list', 'nclusts'};

% savename = [ROOT_DIR, filesep, 'projectnb/busplab/Experiments/ECoG_Preprocessed_AM/topelc_data_to_surf_60pct']; 
%     vars_to_save = {'elc', 'elc_cons', 'elc_vow', 'elc_word', 'elc_anytop', 'topelc', 'top_proportion_electrodes',...
%         'clustlist', 'clusterkey','global_clust_num_list', 'nclusts'};

%%
 %%%%%%%%%%%%%%%% organize data %%%%%%%%%%
n_elcs = height(elc);  
falsecol = false(n_elcs,1);
elc = [elc, table(falsecol, falsecol, falsecol, falsecol, 'VariableNames', {'top_cons_coder','top_vow_coder','top_word_coder','top_cons_or_word'})];                              
elc.feat_set = []; % delete this variable which makes the filesize much larger

% create table with top-30% coding electrodes data by cluster
clustlist = {'PtM-s','PtM-r','ME-sb','ME-sn','AP-r','AP-s'}';
nclusts = length(clustlist); 

% make tables of top coders
top_inds = round(n_elcs*[1-ops.top_proportion_electrodes]) : n_elcs; 
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
elc.anytop(elc.top_cons_coder | elc.top_vow_coder | elc.top_word_coder) = true; 
    elc_anytop= elc(elc.anytop, :); 
elc.top_cons_or_word(elc.top_cons_coder | elc.top_word_coder) = true; 
    elc_cons_or_word = elc(elc.top_cons_or_word,:); 

elc = sortrows(elc,{'subject','electrode'}); % sort electrodes by subject and electrode ID
    
nancol = NaN(nclusts, 1);
nan2 = [nancol, nancol]; 
topelc = table(clustlist, nancol, nan2,   nancol,         nancol,              nan2,       nan2,          nancol,       nancol,            nan2,      nan2,         nancol,   nancol,           nan2,      nan2,        nancol,          nancol,              nancol,          nancol,  'VariableNames',...
                  {'clust',  'n_elc',  'n_elc_RL', 'cons_prop', 'cons_left_prop', 'cons_RL', 'cons_ebar', 'vow_prop','vow_left_prop', 'vow_RL', 'vow_ebar', 'word_prop','word_left_prop', 'word_RL', 'word_ebar', 'anytop_left_n', 'anytop_left_prop', 'cons_word_left_n', 'cons_word_left_prop'}); 
    
% _prop means the proportion of elc in this cluster which are top encoders
% _left_prop means the proportion of top encoders in this cluster which are in left hemisphere
ebar_lims = NaN(nclusts,2);
alpha= 0.00001 : 0.00001 : 0.99999;  % for errorbars 
for iclust = 1:nclusts
    thisclust = clustlist{iclust};
    clustrows = strcmp(elc.cluster_name, thisclust);
    topelc.n_elc(iclust) = nnz(clustrows);
    topelc.n_elc_RL(iclust,1) = nnz(clustrows & ~isnan(elc.right_hemi));  % right this clust
    topelc.n_elc_RL(iclust,2) = nnz(clustrows & isnan(elc.right_hemi));  % left this clust
    
    % consonant 
    clustrows_top_cons = clustrows & elc.top_cons_coder; 
    n_top_cons_this_clust = nnz(clustrows_top_cons); % number of electrodes in this cluster which are top cons code
    n_top_cons_this_clust_leftside = nnz(clustrows_top_cons & isnan(elc.right_hemi)); % number electrodes in this clust which are top coders and in L hem 
    topelc.cons_RL(iclust,:) = [n_top_cons_this_clust-n_top_cons_this_clust_leftside  , n_top_cons_this_clust_leftside]; % number of R and L electrodes
    topelc.cons_left_prop(iclust) = n_top_cons_this_clust_leftside / n_top_cons_this_clust; 
    topelc.cons_prop(iclust) = nnz(clustrows_top_cons) / topelc.n_elc(iclust); 
    
        % error bars; use binomial distribution; equation from Alfonso Nieto-Castanon
        X = sum(topelc.cons_RL(iclust,:)); % number of top elcs from this clust
        p=binocdf(X,topelc.n_elc(iclust),alpha); 
        topelc.cons_ebar(iclust,:) = alpha([find(p>.975,1,'last'),find(p<.025,1,'first')]); 

    % vowel
    clustrows_top_vow = clustrows & elc.top_vow_coder; 
    n_top_vow_this_clust = nnz(clustrows_top_vow); % number of electrodes in this cluster which are top coder
    n_top_vow_this_clust_leftside = nnz(clustrows_top_vow & isnan(elc.right_hemi)); % number electrodes in this clust which are top coders and in L hem 
    topelc.vow_RL(iclust,:) = [n_top_vow_this_clust-n_top_vow_this_clust_leftside  , n_top_vow_this_clust_leftside]; % number of R and L electrodes
    topelc.vow_left_prop(iclust) = n_top_vow_this_clust_leftside / n_top_vow_this_clust; 
    topelc.vow_prop(iclust) = nnz(clustrows_top_vow) / topelc.n_elc(iclust);

            % error bars; use binomial distribution; equation from Alfonso Nieto-Castanon
        X = sum(topelc.vow_RL(iclust,:)); % number of top elcs from this clust
        p=binocdf(X,topelc.n_elc(iclust),alpha); 
        topelc.vow_ebar(iclust,:) = alpha([find(p>.975,1,'last'),find(p<.025,1,'first')]); 
    
    % syllable
    clustrows_top_word = clustrows & elc.top_word_coder; 
    n_top_word_this_clust = nnz(clustrows_top_word); % number of electrodes in this cluster which are top coder
    n_top_word_this_clust_leftside = nnz(clustrows_top_word & isnan(elc.right_hemi));  
    topelc.word_RL(iclust,:) = [n_top_word_this_clust-n_top_word_this_clust_leftside  , n_top_word_this_clust_leftside]; % number of R and L electrodes
    topelc.word_left_prop(iclust) = n_top_word_this_clust_leftside / n_top_word_this_clust; 
    topelc.word_prop(iclust) = nnz(clustrows_top_word) / topelc.n_elc(iclust);

            % error bars; use binomial distribution; equation from Alfonso Nieto-Castanon
        X = sum(topelc.word_RL(iclust,:)); % number of top elcs from this clust
        p=binocdf(X,topelc.n_elc(iclust),alpha); 
        topelc.word_ebar(iclust,:) = alpha([find(p>.975,1,'last'),find(p<.025,1,'first')]); 
 
    % electrodes which top [cons OR vowel OR syllable]
    clustrows_anytop = clustrows & elc.anytop; 
    n_anytop_this_clust = nnz(clustrows_anytop); % number of electrodes in this cluster which are top coder of any kind
    topelc.anytop_n(iclust,1) = n_anytop_this_clust; 
    n_anytop_this_clust_leftside = nnz(clustrows_anytop & isnan(elc.right_hemi));
    topelc.anytop_left_n(iclust) = n_anytop_this_clust_leftside; 
    topelc.anytop_left_prop(iclust) = n_anytop_this_clust_leftside / n_anytop_this_clust; 

    % electrodes which are top [cons OR syllable]
    clustrows_top_cons_or_word = clustrows & elc.top_cons_or_word; 
    n_top_cons_or_word_this_clust = nnz(clustrows_top_cons_or_word); % number of electrodes in this cluster which are top coder
    topelc.cons_word_n(iclust,1) = n_top_cons_or_word_this_clust; 
    n_top_cons_or_word_this_clust_leftside = nnz(clustrows_top_cons_or_word & isnan(elc.right_hemi));
    topelc.cons_word_left_n(iclust) = n_top_cons_or_word_this_clust_leftside; 
    topelc.cons_word_left_prop(iclust) = n_top_cons_or_word_this_clust_leftside / n_top_cons_or_word_this_clust; 
end

%% add clust timing data
topelc.start = nan(height(topelc),1);
topelc.onset = nan(height(topelc),1);
topelc.peak = nan(height(topelc),1);
topelc.offset = nan(height(topelc),1);
topelc.end = nan(height(topelc),1);

vars_to_copy = {'start','onset','peak','offset','end'};
for iclust = 1:nclusts
   clustname = topelc.clust{iclust};
    matchrow =  strcmp(clustname, clusterkey.name) & strcmp('voice_onset', clusterkey.align_name);
    topelc{iclust,vars_to_copy} = clusterkey{matchrow,vars_to_copy}; 
end
topelc.width = topelc.end - topelc.start;
topelc.width_to_peak = topelc.peak - topelc.start;
topelc.width_from_peak = topelc.end - topelc.peak; 

%%
top_proportion_electrodes = ops.top_proportion_electrodes;
save(savename, vars_to_save{:})