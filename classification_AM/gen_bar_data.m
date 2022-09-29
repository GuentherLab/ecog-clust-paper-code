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