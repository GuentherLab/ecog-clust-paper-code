function plot_sig_electrodes_EZmode(threshd_auc_data_with_mni, sub_num, class_label, alignment)

ovr_labels = cellstr(threshd_auc_data_with_mni.ovr_label); 

switch class_label
    case 'word_name'
        class_ovr_idx_range = 1:12;
    case 'consonants_name'
        class_ovr_idx_range = 13:16;
    case 'vowel_name'
        class_ovr_idx_range = 17:19;
end

data_col = sprintf('%s_sig_data',alignment); 
mni_col = sprintf('%s_mni_coords',alignment);
cluster_col = sprintf('%s_cluster',alignment);

alignment_data = table();

for auc_data_idx = class_ovr_idx_range
    alignment_data = [alignment_data; threshd_auc_data_with_mni.(data_col){auc_data_idx}];
end

alignment_mni = alignment_data{alignment_data.sub_num == sub_num & ~(alignment_data.(cluster_col) == 'NaN'), mni_col};
alignment_clusters = alignment_data{alignment_data.sub_num == sub_num & ~(alignment_data.(cluster_col) == 'NaN'), cluster_col};
alignment_clusters_cell = cellstr(alignment_clusters);

[cluster_color_order, ~, cluster_color_idx] = unique(alignment_clusters_cell);

cm_init = jet; 
ops.colormap = cm_init; 
plot_mni_coords(alignment_mni, cluster_color_idx, ops)
hold on;
cm_final = colormap(jet(length(cluster_color_order)));
cb = colorbar();
cb.Ticks = 1/length(cluster_color_order)/2:1/length(cluster_color_order):1;
cb.TickLabels = cluster_color_order;

title(sprintf('Electrode locations and cluster labels for %s: %s, %s data', strcat('S', string(sub_num)), plaintext(class_label), alignment));
hold off; 
% clear ops cm_init cm_final cb
end
