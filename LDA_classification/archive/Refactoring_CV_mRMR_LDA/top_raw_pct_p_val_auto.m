% close all; 
% clear all; 
clc;
set(0,'DefaultFigureWindowStyle','docked')

if ~exist('p', 'var')
    p = parametersClass('/projectnb/busplab/Experiments/ECoG_fMRI_RS/Experiments/ECoG_Preprocessed_LJ/',... % path to data
        2000,...            % event duration (ms)
        50,...              % window (ms)
        25,...              % stride (ms)
        'none',...    % grouping variable
        150,...             % topN feat pooled
        30)                 % topN feat indiv
end
if ~exist('edb', 'var') 
    [~, edb] = p.generate_database; 
end
if ~exist('edb_nonan', 'var') 
    generate_edb_nonan
end

k_pct = 0.20;

all_stim_counts = groupsummary(edb_nonan(edb_nonan.type == 2, :), {'cluster_name'},...
    'IncludeEmptyGroups', true, 'IncludeMissingGroups', true);
all_onset_counts = groupsummary(edb_nonan(edb_nonan.type == 1, :), {'cluster_name'},...
    'IncludeEmptyGroups', true, 'IncludeMissingGroups', true);

class_labels = p.class_names; 
data_cols = {'word_accuracy_change_wo_electrode',...
    'cons_accuracy_change_wo_electrode',...
    'vowel_accuracy_change_wo_electrode'};

%%

% run_this = false;
% if run_this
%     %%%% WORKING ------ NO TOUCH
%     per_sub_flag = true;
%     swc_per_sub = table();
%     swc_per_sub.cluster_name = unique(edb_nonan.cluster_name);
%     scc_per_sub = swc_per_sub;
%     svc_per_sub = swc_per_sub; 
% 
%     owc_per_sub = swc_per_sub; 
%     occ_per_sub = swc_per_sub; 
%     ovc_per_sub = swc_per_sub;
% 
%     for sub_num = p.sub_nums
%         [sw_sorted_per_sub, ow_sorted_per_sub] = sort_by_sub(edb_nonan, sub_num, 'word_name');
% 
%         [sc_sorted_per_sub, oc_sorted_per_sub] = sort_by_sub(edb_nonan, sub_num, 'consonants_name');
% 
%         [sv_sorted_per_sub, ov_sorted_per_sub] = sort_by_sub(edb_nonan, sub_num, 'vowel_name');
% 
%         [sw_top_k_per_sub, ow_top_k_per_sub] = count_top_k(sw_sorted_per_sub, ow_sorted_per_sub, k_pct); 
%         [sc_top_k_per_sub, oc_top_k_per_sub] = count_top_k(sc_sorted_per_sub, oc_sorted_per_sub, k_pct); 
%         [sv_top_k_per_sub, ov_top_k_per_sub] = count_top_k(sv_sorted_per_sub, ov_sorted_per_sub, k_pct); 
% 
%         swc_per_sub = addvars(swc_per_sub, sw_top_k_per_sub.GroupCount, 'NewVariableNames', {sprintf('S%d', sub_num)});
%         scc_per_sub = addvars(scc_per_sub, sc_top_k_per_sub.GroupCount, 'NewVariableNames', {sprintf('S%d', sub_num)});
%         svc_per_sub = addvars(svc_per_sub, sv_top_k_per_sub.GroupCount, 'NewVariableNames', {sprintf('S%d', sub_num)});
% 
%         owc_per_sub = addvars(owc_per_sub, ow_top_k_per_sub.GroupCount, 'NewVariableNames', {sprintf('S%d', sub_num)});
%         occ_per_sub = addvars(occ_per_sub, oc_top_k_per_sub.GroupCount, 'NewVariableNames', {sprintf('S%d', sub_num)});
%         ovc_per_sub = addvars(ovc_per_sub, ov_top_k_per_sub.GroupCount, 'NewVariableNames', {sprintf('S%d', sub_num)});
% 
%     end
% 
%     sowc_per_sub = table({swc_per_sub, owc_per_sub}, 'VariableNames', {'word_name'});
%     socc_per_sub = table({scc_per_sub, occ_per_sub}, 'VariableNames', {'consonants_name'});
%     sovc_per_sub = table({svc_per_sub, ovc_per_sub}, 'VariableNames', {'vowel_name'});
% 
%     [w_fig_per_sub, sw_bar_data_per_sub, ow_bar_data_per_sub] = plot_props(edb_nonan, sowc_per_sub, k_pct, p, per_sub_flag);
%     [c_fig_per_sub, sc_bar_data_per_sub, oc_bar_data_per_sub] = plot_props(edb_nonan, socc_per_sub, k_pct, p, per_sub_flag);
%     [v_fig_per_sub, sv_bar_data_per_sub, ov_bar_data_per_sub] = plot_props(edb_nonan, sovc_per_sub, k_pct, p, per_sub_flag);
% 
%     stim_prop_data_per_sub = multiouterjoin(sw_bar_data_per_sub.prop_data, sc_bar_data_per_sub.prop_data, sv_bar_data_per_sub.prop_data);
%     onset_prop_data_per_sub = multiouterjoin(ow_bar_data_per_sub.prop_data, oc_bar_data_per_sub.prop_data, ov_bar_data_per_sub.prop_data);
% 
%     p_values_per_sub = chitest_p_values(stim_prop_data_per_sub, onset_prop_data_per_sub, k_pct)
% 
% end

%%%%%%%%%
%%
%%%% WORKING ------ NO TOUCH
k_pcts = [0.1 0.2 0.3 0.4];
for k_pct = k_pcts
per_sub_flag = false;
swc = table();
swc.cluster_name = unique(edb_nonan.cluster_name);
scc = swc;
svc = swc;

owc = swc; 
occ = swc; 
ovc = swc;

class_label = 'word_name';
[sw_sorted, ow_sorted] = sort_by_score(edb_nonan, class_label);
% generate_all_surf_cases(sw_sorted, ow_sorted, class_label, k_pct);

class_label = 'consonants_name';
[sc_sorted, oc_sorted] = sort_by_score(edb_nonan, class_label);
% generate_all_surf_cases(sc_sorted, oc_sorted, class_label, k_pct);

class_label = 'vowel_name';
[sv_sorted, ov_sorted] = sort_by_score(edb_nonan, class_label);
% generate_all_surf_cases(sv_sorted, ov_sorted, class_label, k_pct);

[sw_top_k, ow_top_k] = count_top_k(sw_sorted, ow_sorted, k_pct); 
[sc_top_k, oc_top_k] = count_top_k(sc_sorted, oc_sorted, k_pct); 
[sv_top_k, ov_top_k] = count_top_k(sv_sorted, ov_sorted, k_pct);

for sub_num = p.sub_nums
    swc = add_counts(swc, sw_top_k, sub_num);
    scc = add_counts(scc, sc_top_k, sub_num);
    svc = add_counts(svc, sv_top_k, sub_num);
    
    owc = add_counts(owc, ow_top_k, sub_num);
    occ = add_counts(occ, oc_top_k, sub_num);
    ovc = add_counts(ovc, ov_top_k, sub_num);
end

sowc = table({swc, owc}, 'VariableNames', {'word_name'});
socc = table({scc, occ}, 'VariableNames', {'consonants_name'});
sovc = table({svc, ovc}, 'VariableNames', {'vowel_name'});

[w_fig, sw_bar_data, ow_bar_data] = plot_props(edb_nonan, sowc, k_pct, p, per_sub_flag);
[c_fig, sc_bar_data, oc_bar_data] = plot_props(edb_nonan, socc, k_pct, p, per_sub_flag);
[v_fig, sv_bar_data, ov_bar_data] = plot_props(edb_nonan, sovc, k_pct, p, per_sub_flag);

stim_prop_data = multiouterjoin(sw_bar_data.prop_data, sc_bar_data.prop_data, sv_bar_data.prop_data);
onset_prop_data = multiouterjoin(ow_bar_data.prop_data, oc_bar_data.prop_data, ov_bar_data.prop_data);

p_values_not_sub_sorted = chitest_p_values(stim_prop_data, onset_prop_data, k_pct)
end
%  close all; 

function [sub_stim_sorted, sub_onset_sorted] = sort_by_sub(edb_nonan, sub_num, class_label)

data_cols = {'word_accuracy_change_wo_electrode',...
    'cons_accuracy_change_wo_electrode',...
    'vowel_accuracy_change_wo_electrode'};
switch class_label
    case 'word_name'
        data_col = 1; 
    case 'consonants_name'
        data_col = 2; 
    case 'vowel_name'
        data_col = 3; 
end

sub_stim_sorted = sortrows(edb_nonan(edb_nonan.type == 2 & edb_nonan.subject == sub_num,...
    :),...
    {'subject', data_cols{data_col}},...
    {'ascend', 'descend'});
sub_onset_sorted = sortrows(edb_nonan(edb_nonan.type == 1 & edb_nonan.subject == sub_num,...
    :),...
    {'subject', data_cols{data_col}},...
    {'ascend', 'descend'});

end

function [sub_stim_sorted, sub_onset_sorted] = sort_by_score(edb_nonan, class_label)

data_cols = {'word_accuracy_change_wo_electrode',...
    'cons_accuracy_change_wo_electrode',...
    'vowel_accuracy_change_wo_electrode'};
switch class_label
    case 'word_name'
        data_col = 1; 
    case 'consonants_name'
        data_col = 2; 
    case 'vowel_name'
        data_col = 3; 
end

sub_stim_sorted = sortrows(edb_nonan(edb_nonan.type == 2,...
    :),...
    {data_cols{data_col}},...
    {'descend'});
sub_onset_sorted = sortrows(edb_nonan(edb_nonan.type == 1,...
    :),...
    {data_cols{data_col}},...
    {'descend'});

end

function [stim_top_k_counts, onset_top_k_counts] = count_top_k(stim_sorted, onset_sorted, k_pct)

stim_top_k_counts = groupsummary(stim_sorted(1:round(height(stim_sorted) * k_pct), :),...
    {'subject', 'cluster_name'}, 'IncludeEmptyGroups', true, 'IncludeMissingGroups', true);
onset_top_k_counts = groupsummary(onset_sorted(1:round(height(onset_sorted) * k_pct), :),...
    {'subject', 'cluster_name'}, 'IncludeEmptyGroups', true, 'IncludeMissingGroups', true);
end

function combined = add_counts(counts_tally, top_k_counts, sub_num)

count_vector = top_k_counts.GroupCount(top_k_counts.subject == sub_num);
if ~isempty(count_vector)
    combined = addvars(counts_tally, count_vector, 'NewVariableNames', {sprintf('S%d', sub_num)});
else
    count_vector = zeros(size(counts_tally(:,1)));
    combined = addvars(counts_tally, count_vector, 'NewVariableNames', {sprintf('S%d', sub_num)});
end
end

function [clust_prop_fig, stim_bar_data, onset_bar_data] = plot_props(edb_nonan, stim_onset_counts, k_pct, params, per_sub_flag)

title_size = 28;
subtitle_size = 22; 
label_size = 24; 
legend_size = 18; 
tick_size = 18; 

class_label = stim_onset_counts.Properties.VariableNames{1};

switch class_label
    case 'word_name'
        class_str = 'Word Name'; 
    case 'consonants_name'
        class_str = 'Consonants Name'; 
    case 'vowel_name'
        class_str = 'Vowel Name';
end

stim_all_subs_counts_summary = groupsummary(edb_nonan(edb_nonan.type == 2, :), {'cluster_name'},...
    'IncludeEmptyGroups', true, 'IncludeMissingGroups', true);
onset_all_subs_counts_summary = groupsummary(edb_nonan(edb_nonan.type == 1, :), {'cluster_name'},...
    'IncludeEmptyGroups', true, 'IncludeMissingGroups', true);

stim_per_sub_counts = stim_onset_counts.(class_label){1};
onset_per_sub_counts = stim_onset_counts.(class_label){2}; 

stim_rows = [1, 2, 4, 5, 6, 8];
non_stim_rows = [3, 7]; 
stim_nanbars = stim_per_sub_counts.cluster_name(non_stim_rows);

onset_rows = 3:8;
non_onset_rows = [1, 2];
onset_nanbars = onset_per_sub_counts.cluster_name(non_onset_rows);

xmark_offset = 0.1;

clust_prop_fig = figure(); 
if per_sub_flag
    sgtitle(sprintf('Cluster Proportions in Top %2.0f%% of %s Scores, Per Subject', k_pct * 100, class_str), 'FontSize', title_size)
else 
    sgtitle({'Cluster Proportions in ', sprintf('Top %2.0f%% of All', k_pct * 100), sprintf('%s Scores', class_str)}, 'FontSize', title_size)
end

subplot(2,1,1);
stim_bar_data = gen_bar_data(stim_all_subs_counts_summary, stim_per_sub_counts, class_label);
bar(stim_bar_data.clusters, stim_bar_data.per_sub_props, 'stacked');

hold on;

s_counts_top_k = stim_bar_data.prop_data{stim_rows, 2};
s_counts_all = stim_bar_data.prop_data{stim_rows, 3};
s_test = chitest(s_counts_top_k, s_counts_all, k_pct);

yline(k_pct);
scatter([stim_nanbars(:); stim_nanbars(:)],...
    [repmat(k_pct + xmark_offset,length(stim_nanbars),1); 
    repmat(k_pct - xmark_offset,length(stim_nanbars),1)],...
    'x','k', 'LineWidth', 2.5);
ylim([0,0.7]);
title('Stim', 'FontSize', subtitle_size);
ylabel(sprintf('Proportion in top %2.0f%%', k_pct * 100))
set(gca, 'FontSize', 20);
leg = legend(params.sub_nums_formal, 'Location', 'eastoutside', 'FontSize', legend_size);
an = annotation('textbox', 'Position', [leg.Position(1) leg.Position(2)-.1 .4 .4], 'String', sprintf('p-val = %0.3f', s_test), 'FitBoxToText', 'on');
an.FontSize = legend_size;
an.FontWeight = 'bold';
if s_test <= 0.05
    an.BackgroundColor = [0.9290 0.6940 0.1250]
end
hold off; 

subplot(2,1,2);
onset_bar_data = gen_bar_data(onset_all_subs_counts_summary, onset_per_sub_counts, class_label);
bar(onset_bar_data.clusters, onset_bar_data.per_sub_props, 'stacked');

hold on;

o_counts_top_k = onset_bar_data.prop_data{onset_rows, 2};
o_counts_all = onset_bar_data.prop_data{onset_rows, 3};
o_test = chitest(o_counts_top_k, o_counts_all, k_pct);

yline(k_pct);
scatter([onset_nanbars(:); onset_nanbars(:)],...
    [repmat(k_pct + xmark_offset,length(onset_nanbars),1); 
    repmat(k_pct - xmark_offset,length(onset_nanbars),1)],...
    'x','k', 'LineWidth', 2.5);
ylim([0,0.7]);
title('Onset', 'FontSize', subtitle_size);
ylabel(sprintf('Proportion in top %2.0f%%', k_pct * 100))
set(gca, 'FontSize', 20);
leg = legend(params.sub_nums_formal, 'Location', 'eastoutside','FontSize', legend_size);
an = annotation('textbox', 'Position', [leg.Position(1) leg.Position(2)-.1 .4 .4], 'String', sprintf('p-val = %0.3f', o_test), 'FitBoxToText', 'on');
an.FontSize = legend_size;
an.FontWeight = 'bold';
if o_test <= 0.05
    an.BackgroundColor = [0.9290 0.6940 0.1250]
end
hold off;

function bar_data = gen_bar_data(all_subs_counts_summary, per_sub_counts, class_label)

bar_data = struct();
bar_data.clusters = per_sub_counts{:,1};
bar_data.all_subs_counts = all_subs_counts_summary.GroupCount;
bar_data.per_sub_top_k_counts = per_sub_counts{:,2:end}; 
bar_data.all_subs_top_k_counts = sum(bar_data.per_sub_top_k_counts, 2);
bar_data.per_sub_props = bar_data.per_sub_top_k_counts ./ bar_data.all_subs_counts;
bar_data.all_subs_props = sum(bar_data.per_sub_props, 2);

bar_data.prop_data = table(bar_data.clusters,...
    bar_data.all_subs_top_k_counts,...
    bar_data.all_subs_counts,...
    bar_data.all_subs_props, ...
    'VariableNames', {'cluster_name', [class_label, '_counts_top_k'],...
    [class_label, '_counts_all'], [class_label, '_props']});
end

end

function joined = multiouterjoin(varargin)
joined = varargin{1};
    for k = 2:nargin
      joined = outerjoin(joined, varargin{k}, 'MergeKeys', true);
    end
end

function out = chitest_p_values(stim_prop_data, onset_prop_data, k_pct)
top_elcs_proportion = k_pct; % proportion of electrodes which were chosen as the 'best performing' set
onset_rows = 3:8;
stim_rows = [1, 2, 4, 5, 6, 8];

stim = stim_prop_data;
onset = onset_prop_data; 

if ~isempty(stim_prop_data) & isempty(onset_prop_data)
p_stim_word = chitest(stim.word_name_counts_top_k(stim_rows), stim.word_name_counts_all(stim_rows), top_elcs_proportion);
p_stim_consonant = chitest(stim.consonants_name_counts_top_k(stim_rows), stim.consonants_name_counts_all(stim_rows), top_elcs_proportion);
p_stim_vowel = chitest(stim.vowel_name_counts_top_k(stim_rows), stim.vowel_name_counts_all(stim_rows), top_elcs_proportion);

out = table(p_stim_word, p_stim_consonant, p_stim_vowel,...
    'VariableNames', {'p_stim_word', 'p_stim_consonant', 'p_stim_vowel'});
end

if ~isempty(onset_prop_data) & isempty(stim_prop_data)
p_onset_word = chitest(onset.word_name_counts_top_k(onset_rows), onset.word_name_counts_all(onset_rows), top_elcs_proportion);
p_onset_consonant = chitest(onset.consonants_name_counts_top_k(onset_rows), onset.consonants_name_counts_all(onset_rows), top_elcs_proportion);
p_onset_vowel = chitest(onset.vowel_name_counts_top_k(onset_rows), onset.vowel_name_counts_all(onset_rows), top_elcs_proportion);

out = table(p_onset_word, p_onset_consonant, p_onset_vowel, ...
    'VariableNames', {'p_onset_word', 'p_onset_consonant', 'p_onset_vowel'});
end

if ~isempty(stim_prop_data) & ~isempty(onset_prop_data)
p_stim_word = chitest(stim.word_name_counts_top_k(stim_rows), stim.word_name_counts_all(stim_rows), top_elcs_proportion);
p_stim_consonant = chitest(stim.consonants_name_counts_top_k(stim_rows), stim.consonants_name_counts_all(stim_rows), top_elcs_proportion);
p_stim_vowel = chitest(stim.vowel_name_counts_top_k(stim_rows), stim.vowel_name_counts_all(stim_rows), top_elcs_proportion);

p_onset_word = chitest(onset.word_name_counts_top_k(onset_rows), onset.word_name_counts_all(onset_rows), top_elcs_proportion);
p_onset_consonant = chitest(onset.consonants_name_counts_top_k(onset_rows), onset.consonants_name_counts_all(onset_rows), top_elcs_proportion);
p_onset_vowel = chitest(onset.vowel_name_counts_top_k(onset_rows), onset.vowel_name_counts_all(onset_rows), top_elcs_proportion);

out = table(p_onset_word, p_onset_consonant, p_onset_vowel, p_stim_word, p_stim_consonant, p_stim_vowel,...
    'VariableNames', {'p_onset_word', 'p_onset_consonant', 'p_onset_vowel', 'p_stim_word', 'p_stim_consonant', 'p_stim_vowel'});

end
end

function [chisq_pval] = chitest(top_elcs_per_cluster, total_elcs_per_cluster, top_elcs_proportion)

noutcomes = length(top_elcs_per_cluster);
degfr = noutcomes-1;
expected = top_elcs_proportion * total_elcs_per_cluster; % top electrodes if they were proportionately distributed across clusters
devstat = (top_elcs_per_cluster-expected).^2 ./ expected; % calculate deviations from expected for each outcome
chistat = sum(devstat); % chi2 test statistic
chisq_pval = chi2cdf(chistat,degfr,'upper'); % get pval

end

function top_pct_brain_plots(stable, otable, class_label, k_pct, hemi, view)
ops = struct(); 
ops.force_close = 1; 
ops.clusters_as_colors = 1; 
ops.subjects_as_shapes = 1; 
ops.selected_hemi = hemi;
ops.viewpoint = view;
ops.top_k_pct = k_pct;

ops.img_savename = sprintf('top_%2.0f_pct_stim_%s_%s_hemi_%s_view.png', k_pct* 100, class_label, ops.selected_hemi, ops.viewpoint); 
plot_electrodes_on_brain(stable(1:round(height(stable) * k_pct), :), class_label, ops);

ops.img_savename = sprintf('top_%2.0f_pct_onset_%s_%s_hemi_%s_view.png', k_pct* 100, class_label, ops.selected_hemi, ops.viewpoint); 
plot_electrodes_on_brain(otable(1:round(height(otable) * k_pct), :), class_label, ops);

end

function generate_all_surf_cases(stable, otable, class_label, k_pct)

top_pct_brain_plots(stable, otable, class_label, k_pct, 'left', 'left');
top_pct_brain_plots(stable, otable, class_label, k_pct, 'left', 'medial');
top_pct_brain_plots(stable, otable, class_label, k_pct, 'right', 'right');
top_pct_brain_plots(stable, otable, class_label, k_pct, 'right', 'medial');

end
