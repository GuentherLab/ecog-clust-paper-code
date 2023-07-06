%%%% indiv electrodes, proportions of electrodes found in significant data
%%%% compared to total

title_size = 28;
subtitle_size = 22; 
label_size = 24; 
legend_size = 18; 
tick_size = 18; 

electrodes = p.electrodes_table;
if ~strcmp(p.grouping_variable, 'none')
    total_stim_e = height(electrodes(strcmp(electrodes.(p.grouping_variable), p.current_group_value) & electrodes.type == 2, :));
    total_onset_e = height(electrodes(strcmp(electrodes.(p.grouping_variable), p.current_group_value) & electrodes.type == 1, :));
else
    total_stim_e = height(electrodes(electrodes.type == 2, :));
    total_onset_e = height(electrodes(electrodes.type == 1, :));
end

class_stim_e_bool_col = []; 
class_onset_e_bool_col = []; 
 
for class_label_idx = 1:length(class_labels)
    class_label = class_labels{class_label_idx};
    
    switch class_label
        case 'word_name'
            class_ovr_idx_range = 1:12; 
            class_ovr_labels = word_ovr; 
        case 'consonants_name'
            class_ovr_idx_range = 13:16;
            class_ovr_labels = cons_ovr; 
        case 'vowel_name'
            class_ovr_idx_range = 17:19;
            class_ovr_labels = vowel_ovr; 
    end

    stim_e_sig_prop = icdf('Binomial', 1-p.p_value/length(class_ovr_idx_range), total_stim_e, p.p_value) / total_stim_e;
    
    class_stim_e = cellfun(@height, threshd_auc_data.stim_sig_data(class_ovr_idx_range,:));
    class_stim_e_prop = class_stim_e/total_stim_e; 
    
    class_stim_e_bool = class_stim_e_prop >= stim_e_sig_prop; 
    class_stim_e_bool_col = cat(1, class_stim_e_bool_col, class_stim_e_bool);
    
    onset_e_sig_prop = icdf('Binomial', 1-p.p_value/length(class_ovr_idx_range), total_onset_e, p.p_value) / total_onset_e;
    
    class_onset_e = cellfun(@height, threshd_auc_data.onset_sig_data(class_ovr_idx_range,:));
    class_onset_e_prop = class_onset_e/total_onset_e;

    class_onset_e_bool = class_onset_e_prop >= onset_e_sig_prop; 
    class_onset_e_bool_col = cat(1, class_onset_e_bool_col, class_onset_e_bool);
    
    fig_general_name = sprintf('sig_e_ovr_proportions_%s', class_label);
    fig_quickref_filename = [p.times_label, '_', p.current_group_value, '_', fig_general_name, '.fig'];
    png_quickref_filename = [p.times_label, '_', p.current_group_value, '_', fig_general_name, '.png'];
    png_filename = [fig_general_name, '.png'];

    class_stim_onset_prop_fig = figure('FileName', png_filename);
    class_stim_onset_prop_fig.WindowState = 'maximized';
    class_stim_onset_prop_fig.Units = 'normalized';
    
    sgtitle(sprintf('Ratio of electrodes grouped as %s considered significant\nto total electrodes per alignment condition\nClass Label: %s. %s', plaintext(p.current_group_value), plaintext(class_label), plainformat(p.times_label)),'FontSize', title_size);
%     title(sprintf('Class Label: %s. %s', plaintext(class_label), plainformat(p.times_label)),'FontSize', title_size);
    
    subplot(1, 2, 1);
    ax1 = gca; 
    b1 = bar(class_stim_e_prop); 
    ax1.XTickLabel = class_ovr_labels;
    if strcmp(class_label, 'word_name')
        ax1.XTickLabelRotation = 45; 
    end
    ax1.FontSize = label_size;
    title('Stimulus', 'FontSize', subtitle_size);
    ylim([0, .2]);
    yline(stim_e_sig_prop, '--k', 'LineWidth', 2);    
%     ax1.YTick = sort([p_value, 0:0.02:0.20]);
    ylabel('Proportion', 'FontSize', label_size);
    grid on;
    leg1 = legend({'Proportion', sprintf('Significance (%0.4f)', stim_e_sig_prop)}, 'Location', 'northeast', 'FontSize', legend_size);

    subplot(1, 2, 2);
    ax2 = gca; 
    b2 = bar(class_onset_e_prop);
    ax2.XTickLabel = class_ovr_labels; 
    if strcmp(class_label, 'word_name')
        ax2.XTickLabelRotation = 45; 
    end
    ax2.FontSize = label_size;
    title('Onset', 'FontSize', subtitle_size);
    ylim([0, .2]); 
    yline(onset_e_sig_prop, '--k', 'LineWidth', 2);
%     ax2.YTick = sort([p_value, 0:0.02:0.20]);
    ylabel('Proportion', 'FontSize', label_size);
    grid on;
    leg2 = legend({'Proportion', sprintf('Significance (%0.4f)', onset_e_sig_prop)}, 'Location', 'northeast', 'FontSize', legend_size);

    savefig(fullfile(p.figs_quickref_path, fig_quickref_filename));
    exportgraphics(class_stim_onset_prop_fig, fullfile(p.figs_quickref_path, png_quickref_filename))
    exportgraphics(class_stim_onset_prop_fig, fullfile(figures_path, class_stim_onset_prop_fig.FileName))
end

stim_sig_benf_e_auc = table(threshd_auc_data.ovr_label(logical(class_stim_e_bool_col)),...
    threshd_auc_data.stim_sig_data(logical(class_stim_e_bool_col)),...
    'VariableNames', {'ovr_label', 'stim_sig_bool'})

onset_sig_benf_e_auc = table(threshd_auc_data.ovr_label(logical(class_onset_e_bool_col)),...
    threshd_auc_data.onset_sig_data(logical(class_onset_e_bool_col)),...
    'VariableNames', {'ovr_label', 'onset_sig_bool'})

sig_benf_data = outerjoin(stim_sig_benf_e_auc, onset_sig_benf_e_auc, 'MergeKeys', 1)


