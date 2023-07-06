clust_prop_gen(p, edb_nonan(edb_nonan.type==2,:), k_pct)

function [bar_data, prop_data] = clust_prop_gen(p, stim_or_onset_nonan, k_pct)

data_cols = {'word_accuracy_change_wo_electrode', 'cons_accuracy_change_wo_electrode', 'vowel_accuracy_change_wo_electrode'};

w_props_sum = 0;
c_props_sum = 0;
v_props_sum = 0;

for sub_num = p.sub_nums
    sub_num
    sub_nonan = stim_or_onset_nonan(stim_or_onset_nonan.subject == sub_num,:);
    
    k_pct_rows = 1:round(k_pct * height(sub_nonan)); 
    
    word_sorted = sortrows(sub_nonan, data_cols{1}, 'descend');
    word_sorted_top_k_pct = word_sorted(k_pct_rows, :);
    
    cons_sorted = sortrows(sub_nonan, data_cols{2}, 'descend');
    cons_sorted_top_k_pct = cons_sorted(k_pct_rows, :);
    
    vowel_sorted = sortrows(sub_nonan, data_cols{3}, 'descend');
    vowel_sorted_top_k_pct = vowel_sorted(k_pct_rows, :);    
    
    word_info = groupsummary(word_sorted, 'cluster_name', 'IncludeMissingGroups', true, 'IncludeEmptyGroups', true);
    w_counts_all = word_info.GroupCount
    word_top_k_info = groupsummary(word_sorted_top_k_pct, 'cluster_name', 'IncludeMissingGroups', true, 'IncludeEmptyGroups', true);
    w_counts_top_k_pct = word_top_k_info.GroupCount

    w_props = w_counts_top_k_pct ./ w_counts_all
    w_props_sum = w_props_sum + w_props;
    
    cons_info = groupsummary(cons_sorted, 'cluster_name', 'IncludeMissingGroups', true, 'IncludeEmptyGroups', true);
    c_counts_all = cons_info.GroupCount
    cons_top_k_info = groupsummary(cons_sorted_top_k_pct, 'cluster_name', 'IncludeMissingGroups', true, 'IncludeEmptyGroups', true);
    c_counts_top_k_pct = cons_top_k_info.GroupCount

    c_props = c_counts_top_k_pct ./ c_counts_all
    c_props_sum = c_props_sum + c_props;

    vowel_info = groupsummary(vowel_sorted, 'cluster_name', 'IncludeMissingGroups', true, 'IncludeEmptyGroups', true);
    v_counts_all = vowel_info.GroupCount
    vowel_top_k_info = groupsummary(vowel_sorted_top_k_pct, 'cluster_name', 'IncludeMissingGroups', true, 'IncludeEmptyGroups', true);
    v_counts_top_k_pct = vowel_top_k_info.GroupCount

    v_props = v_counts_top_k_pct ./ v_counts_all
    v_props_sum = v_props_sum + v_props;
    
end
w_props_avg = w_props_sum ./ length(p.sub_nums)
c_props_avg = c_props_sum ./ length(p.sub_nums)
v_props_avg = v_props_sum ./ length(p.sub_nums)

    bar_data = table(word_info.cluster_name, w_props_avg, c_props_avg, v_props_avg, ...
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
