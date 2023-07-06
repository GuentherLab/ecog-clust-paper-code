%%%% all electrodes, one-vs-rest labels 

title_size = 28;
subtitle_size = 22; 
label_size = 24; 
legend_size = 18; 
tick_size = 18; 

all_ovr_data = all_class_sub_cv_mRMR_LDA_OvR_data; 

class_labels = fieldnames(all_ovr_data);
sub_nums = [357, 362, 369, 372, 376]'; 

sub_stim_mean_auc_table = table(sub_nums, 'VariableNames', {'sub_num'});
sub_onset_mean_auc_table = table(sub_nums, 'VariableNames', {'sub_num'});

for class_label_idx = 1:numel(class_labels)
    class_label = class_labels{class_label_idx};
    
    stim_class_all_auc = []; 
    onset_class_all_auc = []; 
    
    for sub_num_idx = 1:length(sub_nums)
        sub_num = sub_nums(sub_num_idx); 
        sub_str = strcat('S', string(sub_num));

        sub_stim_class_all_ovr_data = all_ovr_data.(class_label).(sub_str).stim_data;
        sub_onset_class_all_ovr_data = all_ovr_data.(class_label).(sub_str).onset_data;
    
        class_ovr_labels = sub_stim_class_all_ovr_data.Properties.VariableNames; 
        class_ovr_labels(strcmp(class_ovr_labels, 'fold')) = []; 
        
        for class_ovr_label_idx = 1:length(class_ovr_labels)
            class_ovr_label = class_ovr_labels{class_ovr_label_idx};
            
            stim_class_all_auc(sub_num_idx, class_ovr_label_idx) = mean(sub_stim_class_all_ovr_data.(class_ovr_label).stim_auc);
            onset_class_all_auc(sub_num_idx, class_ovr_label_idx) = mean(sub_onset_class_all_ovr_data.(class_ovr_label).onset_auc);
        end
        stim_class_all_auc_table_temp = array2table(stim_class_all_auc, 'VariableNames', class_ovr_labels);
        onset_class_all_auc_table_temp = array2table(onset_class_all_auc, 'VariableNames', class_ovr_labels);
    end
    sub_stim_mean_auc_table = [sub_stim_mean_auc_table, stim_class_all_auc_table_temp];
    sub_onset_mean_auc_table = [sub_onset_mean_auc_table, onset_class_all_auc_table_temp];

end

ovr_labels = sub_stim_mean_auc_table.Properties.VariableNames;

i = 1; j = 1; k = 1; 
for ovr_label_idx = 1:numel(ovr_labels)
    ovr_label = ovr_labels{ovr_label_idx}; 
    
    if length(ovr_label) == 4
        word_ovr{i} = ovr_label; 
        i = i+1; 
    elseif length(ovr_label) == 2
        cons_ovr{j} = ovr_label;
        j = j+1;
    elseif length(ovr_label) == 1
        vowel_ovr{k} = ovr_label;
        k = k+1; 
    end

end

word_ovr_cat = categorical(word_ovr);
cons_ovr_cat = categorical(cons_ovr);
vowel_ovr_cat = categorical(vowel_ovr);

for class_label_idx = 1:length(class_labels)
    class_label = class_labels{class_label_idx};
    switch class_label
        case 'word_name'
            class_label_ovr = word_ovr;
        case 'consonants_name'
            class_label_ovr = cons_ovr;
        case 'vowel_name'
            class_label_ovr = vowel_ovr;
    end

    %%%% class_label Figure
    
    fig_general_name = sprintf('all_ovr_%s', class_label);
    fig_quickref_filename = [p.times_label, '_', p.current_group_value, '_', fig_general_name, '.fig'];
    png_quickref_filename = [p.times_label, '_', p.current_group_value, '_', fig_general_name, '.png'];
    png_filename = [fig_general_name, '.png'];
    
    class_label_auc_fig = figure('FileName', png_filename); 
    class_label_auc_fig.WindowState = 'maximized'; 

    class_label_auc_ax = gca; 

    stim_class_label_auc_arr = table2array(sub_stim_mean_auc_table(:, class_label_ovr));
    stim_class_label_auc_means = mean(stim_class_label_auc_arr, 1)';
    stim_class_label_auc_std = std(stim_class_label_auc_arr, 0, 1)';

    onset_class_label_auc_arr = table2array(sub_onset_mean_auc_table(:, class_label_ovr));
    onset_class_label_auc_means = mean(onset_class_label_auc_arr, 1)';
    onset_class_label_auc_std = std(onset_class_label_auc_arr, 0, 1)';

    class_label_means = [stim_class_label_auc_means, onset_class_label_auc_means]; 
    class_label_stds = [stim_class_label_auc_std, onset_class_label_auc_std]; 
    class_label_auc_bar = bar(class_label_means, 'grouped');
    hold on; 

    [class_label_ngroups, class_label_nbars] = size(class_label_means);
    class_label_err_x = nan(class_label_nbars, class_label_ngroups);
    for i = 1:class_label_nbars
        class_label_err_x(i,:) = class_label_auc_bar(i).XEndPoints; 
    end

    class_label_auc_err_bar = errorbar(class_label_err_x', class_label_means, class_label_stds,...
        'k', 'LineWidth', 1.5, 'linestyle', 'none');
    chance_line = yline(.5, '--k', 'LineWidth', 2);   

    legend([class_label_auc_bar(1), class_label_auc_bar(2), class_label_auc_err_bar(1), chance_line],...
        {'Stimulus', 'Onset', 'Std Dev', 'Chance'}, 'FontSize', legend_size);
    class_label_auc_ax.XTickLabel = class_label_ovr; 
    if strcmp(class_label, 'word_name')
        class_label_auc_ax.XTickLabelRotation = 45;
    end
    set(class_label_auc_ax, 'TickLength', [0, 0]);
    class_label_auc_ax.FontSize = tick_size; 
    ylim([0 1]);
    grid on;

    xlabel('Class Label', 'FontSize', label_size);
    ylabel('AUC', 'FontSize', label_size);
    title(sprintf('Mean AUCs of Classifiers Trained on Individual Electrodes Grouped as %s', plaintext(p.current_group_value)), 'FontSize', title_size);
    subtitle(sprintf('%s, One-versus-Rest, Averaged across subjects. %s', plaintext(class_label), plainformat(p.times_label)), 'FontSize', subtitle_size);
    hold off; 

    figures_path = fullfile(p.grouping_path, p.current_group_value, 'Figures');
    savefig(fullfile(p.figs_quickref_path, fig_quickref_filename));
    exportgraphics(class_label_auc_fig, fullfile(p.figs_quickref_path, png_quickref_filename))
    exportgraphics(class_label_auc_fig, fullfile(figures_path, class_label_auc_fig.FileName))
end

