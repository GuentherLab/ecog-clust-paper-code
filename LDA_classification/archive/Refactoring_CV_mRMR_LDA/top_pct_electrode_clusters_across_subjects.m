% close all; clc
% test(p, edb_nonan, k_pct);
% 
% function [bar_data, prop_data] = test(p, edb_nonan, k_pct)

title_size = 28;
subtitle_size = 22; 
label_size = 24; 
legend_size = 18; 
tick_size = 18; 

cluster_order = {'ESP-s', 'ESP-r', 'PtM-s', 'PtM-r', 'ME-sb', 'ME-sn', 'AP-r', 'AP-s'};
k_pct = 0.2;
sub_nums = p.sub_nums; 
class_labels = p.class_names;

all_stim_prop = 0;
all_onset_prop = 0;
all_stim_counts = 0;
all_onset_counts = 0; 
sub_spec_stim_counts = table(); 

for sub_num = sub_nums
    
    % for each subject

    sub_stim_prop = table();
    sub_onset_prop = table();
    sub_stim_counts = table();
    sub_onset_counts = table(); 
    
    for class_label = class_labels
        % for each class label (word, cons, vowel)

        %%%% pull subject specific nonan electrodes from edb_nonan and take
        %%%% top 20%
        sub_stim_nonan = filter_by_sub_and_type(edb_nonan, sub_num, 2, class_label{1});
        sub_stim_nonan_top_k = sub_stim_nonan(1:round(height(sub_stim_nonan) * k_pct), :);

        %%%% initialize sub_stim_counts table with variable of cluster name
        sub_stim_counts.cluster_name = categorical(categories(sub_stim_nonan.cluster_name));
        sub_spec_stim_counts.cluster_name = categorical(categories(sub_stim_nonan.cluster_name));
        %%%% calculate the counts of clusters for the subject and class
        %%%% then append to sub_stim_counts
        sub_stim_class_counts = calc_counts(sub_stim_nonan, sub_stim_nonan_top_k, class_label{1});
        sub_stim_counts = outerjoin(sub_stim_counts, sub_stim_class_counts, 'MergeKeys', true);
        
        %%%% repeat for onset
        sub_onset_nonan = filter_by_sub_and_type(edb_nonan, sub_num, 1, class_label{1});
        sub_onset_nonan_top_k = sub_onset_nonan(1:round(height(sub_onset_nonan) * k_pct), :);
        
        sub_onset_counts.cluster_name = categorical(categories(sub_onset_nonan.cluster_name));
        sub_onset_class_counts = calc_counts(sub_onset_nonan, sub_onset_nonan_top_k, class_label{1});
        sub_onset_counts = outerjoin(sub_onset_counts, sub_onset_class_counts, 'MergeKeys', true);
        
    end
    
    sub_spec_stim_counts = [sub_spec_stim_counts, table(sub_stim_counts, 'VariableNames', {sprintf('S%d', sub_num)})]
    
    all_stim_counts = all_stim_counts + sub_stim_counts{:, 2:end};
    all_onset_counts = all_onset_counts + sub_onset_counts{:, 2:end};

end

all_stim_counts_temp = array2table(all_stim_counts ./ length(p.sub_nums));
all_stim_counts = [table(sub_stim_counts.cluster_name, 'VariableNames', {'cluster_name'}), all_stim_counts_temp];
all_stim_counts.Properties.VariableNames = sub_stim_counts.Properties.VariableNames;
[stim_bar_data, andrew_stim] = calc_props(all_stim_counts);

all_onset_counts_temp = array2table(all_onset_counts ./ length(p.sub_nums));
all_onset_counts = [table(sub_onset_counts.cluster_name, 'VariableNames', {'cluster_name'}), all_onset_counts_temp];
all_onset_counts.Properties.VariableNames = sub_onset_counts.Properties.VariableNames;
[onset_bar_data, andrew_onset] = calc_props(all_onset_counts);

disp(andrew_stim)
disp(andrew_onset)
out = chitest_p_values(andrew_stim, andrew_onset)

%%% Cluster proportions
for class_idx = 1:length(p.class_names)
    figure('Name', sprintf('Cluster Proportions, %s', p.class_names_formal{class_idx}));

    sgtitle(sprintf('Proportions of Cluster name in top %d%%, %s', k_pct*100, p.class_names_formal{class_idx}), 'FontSize', title_size);
    subplot(2,1,1); 
    bar(stim_bar_data{:, 1}, stim_bar_data{:, class_idx + 1});
    hold on;
    yline(k_pct);
    ylim([0,0.7]);
    title('Stim', 'FontSize', subtitle_size);
    set(gca, 'FontSize', 20);

    subplot(2,1,2); 
    bar(onset_bar_data{:, 1}, onset_bar_data{:, class_idx + 1});
    hold on;
    yline(k_pct);
    ylim([0,0.7]);
    title('Onset', 'FontSize', subtitle_size);
    set(gca, 'FontSize', 20);

end

% end

function sub_edb_nonan_sorted = filter_by_sub_and_type(edb_nonan, sub_num, type, class_label)

switch class_label
    case 'word_name'
        data_col = 1; 
    case 'consonants_name'
        data_col = 2; 
    case 'vowel_name'
        data_col = 3; 
end

data_cols = {'word_accuracy_change_wo_electrode',...
    'cons_accuracy_change_wo_electrode',...
    'vowel_accuracy_change_wo_electrode'};

sub_edb_nonan_sorted = sortrows(edb_nonan(edb_nonan.subject == sub_num & edb_nonan.type == type, {'subject', 'cluster_name', data_cols{data_col}}), data_cols{data_col}, 'descend'); 

end

function [sub_class_prop_counts] = calc_counts(edb_nonan, edb_nonan_top_k, class_label)

sub_nonan_info = groupsummary(edb_nonan, 'cluster_name', 'IncludeMissingGroups', true, 'IncludeEmptyGroups', true);
sub_nonan_counts_all = sub_nonan_info.GroupCount;
sub_nonan_top_k_info = groupsummary(edb_nonan_top_k, 'cluster_name', 'IncludeMissingGroups', true, 'IncludeEmptyGroups', true);
sub_nonan_counts_top_k_pct = sub_nonan_top_k_info.GroupCount;

sub_class_prop_counts = table(sub_nonan_info.cluster_name, sub_nonan_counts_top_k_pct, sub_nonan_counts_all,...
    'VariableNames',...
    {'cluster_name',...
    [class_label, '_top_k_pct'], [class_label, '_counts_all']});

end

function [class_props, class_counts] = calc_props(all_counts)

word_props = all_counts{:,2} ./ all_counts{:,3};%, 'VariableNames', {'word_props'}) 
cons_props = all_counts{:,4} ./ all_counts{:,5};%, 'VariableNames', {'cons_props'}); 
vowel_props = all_counts{:,6} ./ all_counts{:,7};%, 'VariableNames', {'vowel_props'});

class_counts = all_counts;
class_counts = addvars(class_counts, word_props, 'After', 'word_name_counts_all', 'NewVariableNames', 'word_props');
class_counts = addvars(class_counts, cons_props, 'After', 'consonants_name_counts_all', 'NewVariableNames', 'cons_props');
class_counts = addvars(class_counts, vowel_props, 'After', 'vowel_name_counts_all', 'NewVariableNames', 'vowel_props');

class_props = class_counts(:,[1,4,7,10]);

end

function out = chitest_p_values(andrew_stim, andrew_onset)
top_elcs_proportion = 0.2; % proportion of electrodes which were chosen as the 'best performing' set
onset_rows = 3:8;
stim_rows = [1, 2, 4, 5, 6, 8];

onset = andrew_onset; 
stim = andrew_stim;

p_onset_word = chitest(onset.word_name_top_k_pct(onset_rows), onset.word_name_counts_all(onset_rows), top_elcs_proportion);
p_onset_consonant = chitest(onset.consonants_name_top_k_pct(onset_rows), onset.consonants_name_counts_all(onset_rows), top_elcs_proportion);
p_onset_vowel = chitest(onset.vowel_name_top_k_pct(onset_rows), onset.vowel_name_counts_all(onset_rows), top_elcs_proportion);

p_stim_word = chitest(stim.word_name_top_k_pct(stim_rows), stim.word_name_counts_all(stim_rows), top_elcs_proportion);
p_stim_consonant = chitest(stim.consonants_name_top_k_pct(stim_rows), stim.consonants_name_counts_all(stim_rows), top_elcs_proportion);
p_stim_vowel = chitest(stim.vowel_name_top_k_pct(stim_rows), stim.vowel_name_counts_all(stim_rows), top_elcs_proportion);

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





























