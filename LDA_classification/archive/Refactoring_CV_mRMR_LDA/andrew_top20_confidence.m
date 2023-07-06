 % test whether there are significantly different proportions of top-20% electrodes in different cluster
 %
 % updated 2021/11/23

 top_elcs_proportion = 0.2; % proportion of electrodes which were chosen as the 'best performing' set
 onset_rows = 3:8;
 stim_rows = [1, 2, 4, 5, 6, 8];
 datafile = '/project/busplab/software/ecog/scripts_LJ/Refactoring_CV_mRMR_LDA/andrew_prop_data.mat';
 
 load(datafile)
onset = andrew_onset; clear andrew_onset
 stim = andrew_stim; clear andrew_stim
 
 p_onset_word = chitest(onset.word_counts_top_k_pct(onset_rows), onset.word_counts_all(onset_rows), top_elcs_proportion)
 p_onset_consonant = chitest(onset.cons_counts_top_k_pct(onset_rows), onset.cons_counts_all(onset_rows), top_elcs_proportion)
 p_onset_vowel = chitest(onset.vowel_counts_top_k_pct(onset_rows), onset.vowel_counts_all(onset_rows), top_elcs_proportion)
 
 p_stim_word = chitest(stim.word_counts_top_k_pct(stim_rows), stim.word_counts_all(stim_rows), top_elcs_proportion)
 p_stim_consonant = chitest(stim.cons_counts_top_k_pct(stim_rows), stim.cons_counts_all(stim_rows), top_elcs_proportion)
 p_stim_vowel = chitest(stim.vowel_counts_top_k_pct(stim_rows), stim.vowel_counts_all(stim_rows), top_elcs_proportion)

  
 %% 
 function [chisq_pval] = chitest(top_elcs_per_cluster, total_elcs_per_cluster, top_elcs_proportion)
 
noutcomes = length(top_elcs_per_cluster);
degfr = noutcomes-1;
expected = top_elcs_proportion * total_elcs_per_cluster; % top electrodes if they were proportionately distributed across clusters
devstat = (top_elcs_per_cluster-expected).^2 ./ expected; % calculate deviations from expected for each outcome
chistat = sum(devstat); % chi2 test statistic
chisq_pval = chi2cdf(chistat,degfr,'upper'); % get pval
 
 end