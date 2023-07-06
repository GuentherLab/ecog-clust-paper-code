%%%% significant electrodes only, one-vs-rest labels

title_size = 28;
subtitle_size = 22; 
label_size = 24; 
legend_size = 18; 
tick_size = 18; 

all_class_stim_ovr_auc_mean = table();
all_class_stim_ovr_auc_std = table();

all_class_onset_ovr_auc_mean = table();
all_class_onset_ovr_auc_std = table();

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

    class_stim_ovr_data = threshd_auc_data(class_ovr_idx_range, 'stim_sig_data');
    class_onset_ovr_data = threshd_auc_data(class_ovr_idx_range, 'onset_sig_data');
    
    class_stim_ovr_auc_mean = [];
    class_stim_ovr_auc_std = [];
    
    class_onset_ovr_auc_mean = [];
    class_onset_ovr_auc_std = [];
    
    for ovr_label_idx = 1:length(class_ovr_labels)
        ovr_label = class_ovr_labels{ovr_label_idx};
        
        class_stim_ovr_auc_mean(ovr_label_idx) = mean(class_stim_ovr_data.stim_sig_data{ovr_label_idx}.mean_auc);
        class_stim_ovr_auc_std(ovr_label_idx) = std(class_stim_ovr_data.stim_sig_data{ovr_label_idx}.mean_auc);
        
        class_onset_ovr_auc_mean(ovr_label_idx) = mean(class_onset_ovr_data.onset_sig_data{ovr_label_idx}.mean_auc);
        class_onset_ovr_auc_std(ovr_label_idx) = std(class_onset_ovr_data.onset_sig_data{ovr_label_idx}.mean_auc);
        
    end
    
    class_stim_ovr_auc_mean_table = array2table(class_stim_ovr_auc_mean, 'VariableNames', class_ovr_labels);
    class_stim_ovr_auc_std_table = array2table(class_stim_ovr_auc_std, 'VariableNames', class_ovr_labels);    
    
    all_class_stim_ovr_auc_mean = [all_class_stim_ovr_auc_mean, class_stim_ovr_auc_mean_table];
    all_class_stim_ovr_auc_std = [all_class_stim_ovr_auc_std, class_stim_ovr_auc_std_table];
    
    class_onset_ovr_auc_mean_table = array2table(class_onset_ovr_auc_mean, 'VariableNames', class_ovr_labels);
    class_onset_ovr_auc_std_table = array2table(class_onset_ovr_auc_std, 'VariableNames', class_ovr_labels);    
    
    all_class_onset_ovr_auc_mean = [all_class_onset_ovr_auc_mean, class_onset_ovr_auc_mean_table];
    all_class_onset_ovr_auc_std = [all_class_onset_ovr_auc_std, class_onset_ovr_auc_std_table];    
    
end

all_class_ovr_auc_mean = [all_class_stim_ovr_auc_mean; all_class_onset_ovr_auc_mean];
all_class_ovr_auc_mean.Properties.RowNames = {'stim', 'onset'};

all_class_ovr_auc_std = [all_class_stim_ovr_auc_std; all_class_onset_ovr_auc_std];
all_class_ovr_auc_std.Properties.RowNames = {'stim', 'onset'};

all_class_ovr_auc_mean_arr = table2array(all_class_ovr_auc_mean)';
all_class_ovr_auc_std_arr = table2array(all_class_ovr_auc_std)';

%%%% Admittedly, what follows is pretty sloppy. Copied and pasted from
%%%% all_ovr_visuals.m and adapted from there. A software engineer would
%%%% look at this and laugh

ovr_labels = all_class_ovr_auc_mean.Properties.VariableNames;

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
    fig_general_name = sprintf('sig_e_ovr_%s', class_label);
    fig_quickref_filename = [p.times_label, '_', p.current_group_value, '_', fig_general_name, '.fig'];
    png_quickref_filename = [p.times_label, '_', p.current_group_value, '_', fig_general_name, '.png'];
    png_filename = [fig_general_name, '.png'];

    sig_e_class_label_auc_fig = figure('FileName', png_filename);
    sig_e_class_label_auc_fig.WindowState = 'maximized'; 

    sig_e_class_label_auc_ax = gca; 

    stim_class_label_auc_means = table2array(all_class_ovr_auc_mean('stim', class_label_ovr))';
    stim_class_label_auc_std = table2array(all_class_ovr_auc_std('stim', class_label_ovr))';

    onset_class_label_auc_means = table2array(all_class_ovr_auc_mean('onset', class_label_ovr))';
    onset_class_label_auc_std = table2array(all_class_ovr_auc_std('onset', class_label_ovr))';

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

    legend([class_label_auc_bar(1), class_label_auc_bar(2), class_label_auc_err_bar(1), chance_line], {'Stimulus', 'Onset', 'STD', 'Chance'}, 'FontSize', legend_size);
    sig_e_class_label_auc_ax.XTickLabel = class_label_ovr; 
    if strcmp(class_label, 'word_name')
        sig_e_class_label_auc_ax.XTickLabelRotation = 45;
    end
    set(sig_e_class_label_auc_ax, 'TickLength', [0, 0]);
    sig_e_class_label_auc_ax.FontSize = tick_size; 
    ylim([0 1]);
    grid on;

    xlabel('Class Label', 'FontSize', label_size);
    ylabel('AUC', 'FontSize', label_size);
    title(sprintf('Mean AUCs of Classifiers Corresponding to\nStatistically Significant Electrodes Grouped as %s', p.current_group_value), 'FontSize', title_size);
    subtitle(sprintf('%s, One-versus-Rest, Averaged across subjects. %s', plaintext(class_label), plainformat(p.times_label)), 'FontSize', subtitle_size);
    hold off;

    savefig(fullfile(p.figs_quickref_path, fig_quickref_filename));
    exportgraphics(sig_e_class_label_auc_fig, fullfile(p.figs_quickref_path, png_quickref_filename));
    exportgraphics(sig_e_class_label_auc_fig, fullfile(figures_path, sig_e_class_label_auc_fig.FileName));
end








