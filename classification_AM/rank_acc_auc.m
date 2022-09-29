 %%% rank electrodes by acc and auc amonst all electrodes within their cluster
 %
 % Andrew Meier 2022/2/2
 
 auc_acc_data_file = '/usr2/postdoc/amsmeier/electrodes_acc_auc.mat';
 
 load(auc_acc_data_file); 

vars_to_rank = {'acc_consonants', 'acc_vowel', 'acc_word', 'auc_consonants', 'auc_vowel', 'auc_word'};
nvars_to_rank = length(vars_to_rank);

for ivar = 1:nvars_to_rank
    thisvar = vars_to_rank{ivar};
    thisvar_rank = ['clustrank_', thisvar];
    electrodes.(thisvar_rank) = NaN(height(electrodes),1);
    
    for this_align = unique(electrodes.type)'
        elc_this_align = electrodes.type==this_align;
        unq_glob_clusts = unique(electrodes.global_clust_num(elc_this_align)); % all clusters in this alignment
        nclusts = length(unq_glob_clusts);
        for iclust = 1:nclusts
            this_glob_clust = unq_glob_clusts(iclust);
            elc_this_clust = elc_this_align & electrodes.global_clust_num==this_glob_clust;
            thisclust_rows = find(elc_this_clust); 
            n_elc_this_clust = nnz(elc_this_clust);
            thisclust_table = electrodes(elc_this_clust,:);
            [~, clustrank] = sort(thisclust_table.(thisvar), 'ascend'); % higher rank = higher performance
            rankfrac = [clustrank-1] / [n_elc_this_clust]; % ranks spread from zero to 1/n_elc_this_clust
            electrodes.(thisvar_rank)(thisclust_rows) = rankfrac;
        end
    end
    
end

clearvars('-except','electrodes')

