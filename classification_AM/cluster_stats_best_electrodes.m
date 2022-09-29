 % test whether there are significantly different proportions of top-20% electrodes in different cluster
 %
 % updated 2021/12/03 by AM
 
 clear
% close all

 elc_pooling = 'pooled'; % first pool electrodes from all subjects, then get top X percent electrodes and count cluster labels
%     elc_pooling = 'sub_specific'; % first get top X percent electrodes from each subject, then pool electrodes and count cluster labels
 
 top_elcs_proportion = 0.2; % proportion of electrodes which were chosen as the 'best performing' set
 onset_rows = 3:8;
    non_onset_rows = [1, 2];
 stim_rows = [1, 2, 4, 5, 6, 8];
    non_stim_rows = [3, 7]; 
 n_subjects = 5; 
 cluster_order =  {'ESP-r', 'ESP-s', 'PtM-s', 'PtM-r', 'ME-sb', 'ME-sn', 'AP-r', 'AP-s' }; %  clusters by approx peak time

 datafile_top_elcs_pooled = '/project/busplab/software/ecog/scripts_LJ/Refactoring_CV_mRMR_LDA/andrew_prop_data.mat';
 datafiles_top_elcs_sub_specific = {'/project/busplab/software/ecog/scripts_LJ/Refactoring_CV_mRMR_LDA/andrew_onset.mat'; ...
     '/project/busplab/software/ecog/scripts_LJ/Refactoring_CV_mRMR_LDA/andrew_stim.mat'};
 
switch elc_pooling
    case 'pooled'
         load(datafile_top_elcs_pooled)
         supertitle = 'best electrodes (ignoring subject labels)'; 
    case 'sub_specific'
        load(datafiles_top_elcs_sub_specific{1}, 'andrew_onset')
        supertitle = 'best electrodes (equal % from each subject)'; 
        counts_columns = contains(andrew_onset.Properties.VariableNames, '_name_'); % columns that need to be converted)
        andrew_onset{:,counts_columns} = n_subjects * andrew_onset{:,counts_columns} ; % convert mean across subjects to totals
        load(datafiles_top_elcs_sub_specific{2}, 'andrew_stim')
        counts_columns = contains(andrew_stim.Properties.VariableNames, '_name_'); % columns that need to be converted)
        andrew_stim{:,counts_columns} = n_subjects * andrew_stim{:,counts_columns} ; % convert mean across subjects to totals
end
onset = andrew_onset; clear andrew_onset
stim = andrew_stim; clear andrew_stim
onset.Properties.VariableNames(2:end) = strrep(onset.Properties.VariableNames(2:end),'name_counts','counts') ; % rename variables
onset.Properties.VariableNames(2:end) = strrep(onset.Properties.VariableNames(2:end),'name','counts') ; % rename variables
onset.Properties.VariableNames(2:end) = strrep(onset.Properties.VariableNames(2:end),'consonants','cons') ; % rename variables
onset.Properties.VariableNames(2:end) = strrep(onset.Properties.VariableNames(2:end),'props','proportions') ; % rename variables
 [~, neworder] = ismember(cluster_order, onset.cluster_name);
    onset = onset(neworder, :); % reorder clusters by approx peak time
stim.Properties.VariableNames(2:end) = strrep(stim.Properties.VariableNames(2:end),'name_counts','counts') ; % rename variables
stim.Properties.VariableNames(2:end) = strrep(stim.Properties.VariableNames(2:end),'name','counts') ; % rename variables
stim.Properties.VariableNames(2:end) = strrep(stim.Properties.VariableNames(2:end),'consonants','cons') ; % rename variables
stim.Properties.VariableNames(2:end) = strrep(stim.Properties.VariableNames(2:end),'props','proportions') ; % rename variables
 [~, neworder] = ismember(cluster_order, stim.cluster_name);
    stim = stim(neworder, :); % reorder clusters by approx peak time

 p_onset_word = chitest(onset.word_counts_top_k_pct(onset_rows), onset.word_counts_all(onset_rows), top_elcs_proportion)
 p_onset_consonant = chitest(onset.cons_counts_top_k_pct(onset_rows), onset.cons_counts_all(onset_rows), top_elcs_proportion)
 p_onset_vowel = chitest(onset.vowel_counts_top_k_pct(onset_rows), onset.vowel_counts_all(onset_rows), top_elcs_proportion)
 
 p_stim_word = chitest(stim.word_counts_top_k_pct(stim_rows), stim.word_counts_all(stim_rows), top_elcs_proportion)
 p_stim_consonant = chitest(stim.cons_counts_top_k_pct(stim_rows), stim.cons_counts_all(stim_rows), top_elcs_proportion)
 p_stim_vowel = chitest(stim.vowel_counts_top_k_pct(stim_rows), stim.vowel_counts_all(stim_rows), top_elcs_proportion)

 %% plotting
 figure
hbar = bargr(onset.cons_proportions, 'consonant onset', non_onset_rows, 1, cluster_order, top_elcs_proportion);
hbar = bargr(onset.vowel_proportions, 'vowel onset', non_onset_rows, 2, cluster_order, top_elcs_proportion);
hbar = bargr(onset.word_proportions, 'word onset', non_onset_rows, 3, cluster_order, top_elcs_proportion);
hbar = bargr(stim.cons_proportions, 'consonant stim', non_stim_rows, 4, cluster_order, top_elcs_proportion);
hbar = bargr(stim.vowel_proportions, 'vowel stim', non_stim_rows, 5, cluster_order, top_elcs_proportion);
hbar = bargr(stim.word_proportions, 'word stim', non_stim_rows, 6, cluster_order, top_elcs_proportion);

suptitle(supertitle)

 %% 
 function [chisq_pval] = chitest(top_elcs_per_cluster, total_elcs_per_cluster, top_elcs_proportion)
 
noutcomes = length(top_elcs_per_cluster);
degfr = noutcomes-1;
expected = top_elcs_proportion * total_elcs_per_cluster; % top electrodes if they were proportionately distributed across clusters
devstat = (top_elcs_per_cluster-expected).^2 ./ expected; % calculate deviations from expected for each outcome
chistat = sum(devstat); % chi2 test statistic
chisq_pval = chi2cdf(chistat,degfr,'upper'); % get pval
 
 end
 
 %% 
 function hbar = bargr(yvals, plot_title, nanbars, subplotind, cluster_order, top_elcs_proportion)
 xmark_offset = 0.1; 
 ylimits = [0 0.7]; 
 x_tick_angle = 75; 
 
 subplot(2,3,subplotind)
 hbar = bar(yvals);
 set(gca,'XTickLabels',cluster_order)
 xtickangle(x_tick_angle)
 hold on
 hscat = scatter([nanbars(:); nanbars(:)], ...
     [repmat(top_elcs_proportion + xmark_offset,length(nanbars),1); repmat(top_elcs_proportion - xmark_offset,length(nanbars),1)], 'x','k');
 yline(top_elcs_proportion)
 hold off
 title(plot_title)
 ylim(ylimits)
 
 
 
 end
 