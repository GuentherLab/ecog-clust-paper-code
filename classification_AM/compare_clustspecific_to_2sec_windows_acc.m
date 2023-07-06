 %%%% compare subject-level accuracies for cons, vow, syl when using cluster-specific vs. fixed analysis windows

% also test whether, in 2sec analysis, subject-level decoding across all subjects was significantly difference from chance
% note: this script uses (incorrect) non-emprical chance levels which assume infinite sample size....
% ..... better assessment of chance levels, run chance_level_accuraacy_computations.m.....
% ..... which uses use binomial formula for computing chance levels from combrisson and jerbi 2015

clear
set_paths()

phonunits =             {'cons','vow','syl'};
n_labels_per_phonunit = [4,     3,     12 ];  

dat_2sec = load([DIR_LDA_2SEC_WINDOW filesep 'single_elc_included_trialwise_accuracy.mat']); 
dat_clustspec = load([DIR_LDA_CLUSTER_SPEC_WINDOW filesep 'single_elc_included_trialwise_accuracy.mat']);

clc

% 
chance_level_accuracy_computations()

for iphonunit = 1:length(phonunits)
    thisphonunit = phonunits{iphonunit}
    [~, p_clustspec_vs_2sec] = ttest(dat_2sec.subdata{:,[thisphonunit,'_acc_overall']}, dat_clustspec.subdata{:,[thisphonunit,'_acc_overall']})
    mean_2sec = mean(dat_2sec.subdata{:,[thisphonunit,'_acc_overall']})
    mean_clustspec = mean(dat_clustspec.subdata{:,[thisphonunit,'_acc_overall']});

    [~,p_decoding_above_chance] = ttest(dat_2sec.subdata{:,[thisphonunit,'_acc_overall']} - 1/n_labels_per_phonunit(iphonunit))
end
