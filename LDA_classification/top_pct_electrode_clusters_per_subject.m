% top k percent of electrodes (across subject pool per subject) analysis
% for cluster prevalence
%
% needs edb_cat and edb_nonan

if ~exist('edb_cat', 'var')
    generate_edb_nonan
end

close all;

title_size = 28;
subtitle_size = 22; 
label_size = 24; 
legend_size = 18; 
tick_size = 18; 

cluster_order = {'ESP-s', 'ESP-r', 'PtM-s', 'PtM-r', 'ME-sb', 'ME-sn', 'AP-r', 'AP-s'};
k_pct = 0.2;

if ~exist('edb_nonan', 'var')
    edb_nonan = edb_cat(~(edb_cat.cluster_name == 'NaN'), :);
end

edb_nonan.cluster_name = removecats(edb_nonan.cluster_name, 'NaN');
edb_nonan.cluster_name = categorical(edb_nonan.cluster_name);
edb_nonan.cluster_name = reordercats(edb_nonan.cluster_name, cluster_order);

all_sub_p_vals = table(); 

for sub_num_idx = 1:length(p.sub_nums)
    sub_num = p.sub_nums(sub_num_idx)
    sub_str = p.sub_nums_formal{sub_num_idx};

    stim_nonan = edb_nonan(edb_nonan.type == 2 & edb_nonan.subject == sub_num, :); 
    onset_nonan = edb_nonan(edb_nonan.type == 1 & edb_nonan.subject == sub_num, :); 

    [stim_nonan_bar_data, andrew_stim] = clust_prop_gen(stim_nonan, k_pct);
    stim_nonan_sum = stim_nonan_bar_data{:, 2:4}
    
    [onset_nonan_bar_data, andrew_onset] = clust_prop_gen(onset_nonan, k_pct);
    onset_nonan_bar_data;
    
    sub_p_vals = andrew_calc(andrew_stim, andrew_onset);
    sub_p_vals = [table(sub_num, 'VariableNames', {'sub_num'}), sub_p_vals];
    all_sub_p_vals = [all_sub_p_vals; sub_p_vals];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% sub specific accuracy distributions
    % for sub_num_idx = 1:length(p.sub_nums)
    %     sub_num = p.sub_nums(sub_num_idx);
    %     sub_str = p.sub_nums_formal{sub_num_idx};
    %     
    %     dev_dist_plot_s = plot_deviation_dist(sub_num, stim_nonan);
    %     dev_dist_plot_o = plot_deviation_dist(sub_num, onset_nonan);    
    %     
    % end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Cluster proportions
    for class_idx = 1:length(p.class_names)
        figure('Name', sprintf('%s, Cluster Proportions, %s', sub_str, p.class_names_formal{class_idx}));

        sgtitle(sprintf('%s, Proportions of Cluster name in top %d%%, %s', sub_str, k_pct*100, p.class_names_formal{class_idx}), 'FontSize', title_size);
        subplot(2,1,1); 
        bar(stim_nonan_bar_data{:, 1}, stim_nonan_bar_data{:, class_idx + 1});
        hold on;
        yline(k_pct);
        ylim([0,0.7]);
        title('Stim', 'FontSize', subtitle_size);
        set(gca, 'FontSize', 20);

        subplot(2,1,2); 
        bar(onset_nonan_bar_data{:, 1}, onset_nonan_bar_data{:, class_idx + 1});
        hold on;
        yline(k_pct);
        ylim([0,0.7]);
        title('Onset', 'FontSize', subtitle_size);
        set(gca, 'FontSize', 20);

    end
end

all_sub_p_vals

close all; 

%%
function [bar_data, prop_data] = clust_prop_gen(edb_nonan_subset, k_pct)
k_pct = round(k_pct * height(edb_nonan_subset));

word_sorted = sortrows(edb_nonan_subset, 'word_accuracy_change_wo_electrode', 'descend'); 
word_sorted_top_k_pct = word_sorted(1:k_pct, :);

cons_sorted = sortrows(edb_nonan_subset, 'cons_accuracy_change_wo_electrode', 'descend'); 
cons_sorted_top_k_pct = cons_sorted(1:k_pct, :);

vowel_sorted = sortrows(edb_nonan_subset, 'vowel_accuracy_change_wo_electrode', 'descend'); 
vowel_sorted_top_k_pct = vowel_sorted(1:k_pct, :);

word_info = groupsummary(word_sorted, 'cluster_name', 'IncludeMissingGroups', true, 'IncludeEmptyGroups', true);
w_counts_all = word_info.GroupCount;
word_top_k_info = groupsummary(word_sorted_top_k_pct, 'cluster_name', 'IncludeMissingGroups', true, 'IncludeEmptyGroups', true);
w_counts_top_k_pct = word_top_k_info.GroupCount;

w_props = w_counts_top_k_pct ./ w_counts_all;

cons_info = groupsummary(cons_sorted, 'cluster_name', 'IncludeMissingGroups', true, 'IncludeEmptyGroups', true);
c_counts_all = cons_info.GroupCount;
cons_top_k_info = groupsummary(cons_sorted_top_k_pct, 'cluster_name', 'IncludeMissingGroups', true, 'IncludeEmptyGroups', true);
c_counts_top_k_pct = cons_top_k_info.GroupCount;

c_props = c_counts_top_k_pct ./ c_counts_all;

vowel_info = groupsummary(vowel_sorted, 'cluster_name', 'IncludeMissingGroups', true, 'IncludeEmptyGroups', true);
v_counts_all = vowel_info.GroupCount;
vowel_top_k_info = groupsummary(vowel_sorted_top_k_pct, 'cluster_name', 'IncludeMissingGroups', true, 'IncludeEmptyGroups', true);
v_counts_top_k_pct = vowel_top_k_info.GroupCount;

v_props = v_counts_top_k_pct ./ v_counts_all;

bar_data = table(word_info.cluster_name, w_props, c_props, v_props, ...
    'VariableNames', {'cluster_name', 'word_proportions', 'cons_proportions', 'vowel_proportions'});

prop_data = table(word_info.cluster_name, w_counts_top_k_pct, w_counts_all, w_props,...
    c_counts_top_k_pct, c_counts_all, c_props,...
    v_counts_top_k_pct, v_counts_all, v_props, ...
    'VariableNames',...
    {'cluster_name',...
    'word_counts_top_k_pct', 'word_counts_all', 'word_proportions',...
    'cons_counts_top_k_pct', 'cons_counts_all', 'cons_proportions',...
    'vowel_counts_top_k_pct', 'vowel_counts_all', 'vowel_proportions'});

end

function dev_dist_plot = plot_deviation_dist(sub_num, edb_nonan)

alignment = edb_nonan.type(1);
switch alignment
    case 2
        alignment = 'stim';
    case 1
        alignment = 'onset';
end

title_size = 28;
subtitle_size = 22; 
label_size = 24; 
legend_size = 18; 
tick_size = 18; 

elect_list = edb_nonan.electrode; 
dev_dist_plot = figure('Name', ['S' + string(sub_num) + ' ' + alignment]);

edges = linspace(-.25, .25, 21); 

sgtitle(sprintf('%s, %s, number of electrodes per accuracy deviation range', ['S' + string(sub_num)], alignment), 'FontSize', title_size);

subplot(3,1,1);
histogram(edb_nonan.word_accuracy_change_wo_electrode(edb_nonan.subject == sub_num), 'BinLimits', [-.25, .25], 'BinEdges', edges)
title('Word', 'FontSize', subtitle_size);
set(gca, 'FontSize', 20);

subplot(3,1,2);
histogram(edb_nonan.cons_accuracy_change_wo_electrode(edb_nonan.subject == sub_num), 'BinLimits', [-.25, .25], 'BinEdges', edges)
title('Consonant', 'FontSize', subtitle_size);
ylabel('Count', 'FontSize', label_size);
set(gca, 'FontSize', 20);

subplot(3,1,3);
histogram(edb_nonan.vowel_accuracy_change_wo_electrode(edb_nonan.subject == sub_num), 'BinLimits', [-.25, .25], 'BinEdges', edges)
title('Vowel', 'FontSize', subtitle_size);
xlabel('\Delta Accuracy', 'FontSize', label_size);
set(gca, 'FontSize', 20);

end

function out = andrew_calc(andrew_stim, andrew_onset)
top_elcs_proportion = 0.2; % proportion of electrodes which were chosen as the 'best performing' set
onset_rows = 3:8;
stim_rows = [1, 2, 4, 5, 6, 8];

onset = andrew_onset; 
stim = andrew_stim;

p_onset_word = chitest(onset.word_counts_top_k_pct(onset_rows), onset.word_counts_all(onset_rows), top_elcs_proportion);
p_onset_consonant = chitest(onset.cons_counts_top_k_pct(onset_rows), onset.cons_counts_all(onset_rows), top_elcs_proportion);
p_onset_vowel = chitest(onset.vowel_counts_top_k_pct(onset_rows), onset.vowel_counts_all(onset_rows), top_elcs_proportion);

p_stim_word = chitest(stim.word_counts_top_k_pct(stim_rows), stim.word_counts_all(stim_rows), top_elcs_proportion);
p_stim_consonant = chitest(stim.cons_counts_top_k_pct(stim_rows), stim.cons_counts_all(stim_rows), top_elcs_proportion);
p_stim_vowel = chitest(stim.vowel_counts_top_k_pct(stim_rows), stim.vowel_counts_all(stim_rows), top_elcs_proportion);

out = table(p_onset_word, p_onset_consonant, p_onset_vowel, p_stim_word, p_stim_consonant, p_stim_vowel,...
    'VariableNames', {'p_onset_word', 'p_onset_consonant', 'p_onset_vowel', 'p_stim_word', 'p_stim_consonant', 'p_stim_vowel'});


function [chisq_pval] = chitest(top_elcs_per_cluster, total_elcs_per_cluster, top_elcs_proportion)
noutcomes = length(top_elcs_per_cluster);
degfr = noutcomes-1;
expected = top_elcs_proportion * total_elcs_per_cluster; % top electrodes if they were proportionately distributed across clusters
devstat = (top_elcs_per_cluster-expected).^2 ./ expected; % calculate deviations from expected for each outcome
chistat = sum(devstat); % chi2 test statistic
chisq_pval = chi2cdf(chistat,degfr,'upper'); % get pval

end

end
