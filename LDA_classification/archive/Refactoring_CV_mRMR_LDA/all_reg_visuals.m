%%%% all electrodes, regular labels

title_size = 28;
subtitle_size = 22; 
label_size = 24; 
legend_size = 18; 
tick_size = 18; 
chance_line_size = 2;

all_reg_data = all_class_sub_cv_mRMR_LDA_data; 

class_labels = fieldnames(all_reg_data);
sub_nums = [357, 362, 369, 372, 376]'; 
sub_strs = cell([1, length(sub_nums)]);

sub_stim_mean_acc_table = table(sub_nums, 'VariableNames', {'sub_num'});
sub_onset_mean_acc_table = table(sub_nums, 'VariableNames', {'sub_num'});

for class_label_idx = 1:numel(class_labels)
    class_label = class_labels{class_label_idx};
    
    sub_stim_mean_acc_table_temp = table(); 
    sub_onset_mean_acc_table_temp = table();
    
    stim_acc_arr_temp = zeros([length(sub_nums), 1]); 
    onset_acc_arr_temp = zeros([length(sub_nums), 1]); 
    for sub_num_idx = 1:length(sub_nums)
        sub_num = sub_nums(sub_num_idx); 
        sub_str = strcat('S', string(sub_num));
        sub_strs{sub_num_idx} = convertStringsToChars(sub_str);

        sub_stim_class_data = all_reg_data.(class_label).(sub_str).stim_data.(class_label);
        sub_onset_class_data = all_reg_data.(class_label).(sub_str).onset_data.(class_label);
    
        stim_acc_arr_temp(sub_num_idx) = mean(sub_stim_class_data.stim_acc);
        onset_acc_arr_temp(sub_num_idx) = mean(sub_onset_class_data.onset_acc);
        
        if sub_num_idx == length(sub_nums)
            stim_acc_table_temp = table(sub_nums, stim_acc_arr_temp, 'VariableNames', {'sub_num', class_label}); 
            onset_acc_table_temp = table(sub_nums, onset_acc_arr_temp, 'VariableNames', {'sub_num', class_label}); 
            
            sub_stim_mean_acc_table = outerjoin(sub_stim_mean_acc_table, stim_acc_table_temp, 'MergeKeys', 1);
            sub_onset_mean_acc_table = outerjoin(sub_onset_mean_acc_table, onset_acc_table_temp, 'MergeKeys', 1); 
        end
    end
end

subject_labels = categorical(sub_strs);

fig_general_name = ['inter_sub_acc'];
fig_quickref_filename = [p.times_label, '_', p.current_group_value, '_', fig_general_name, '.fig'];
png_quickref_filename = [p.times_label, '_', p.current_group_value, '_', fig_general_name, '.png'];
png_filename = [fig_general_name, '.png'];

inter_sub_acc_fig = figure('FileName', png_filename);
inter_sub_acc_fig.WindowState = 'maximized';

sgtitle(sprintf('Inter-subject Comparison of LDA Model Accuracies\nWhen Trained on Multiclass Data, Electrodes Grouped as %s. %s', plaintext(p.current_group_value), plainformat(p.times_label)), 'FontSize', title_size);

inter_sub_acc_stim_ax = subplot(1, 2, 1);
hold on;
stim_acc_arr = 100.*table2array(sub_stim_mean_acc_table(:,{'word_name', 'consonants_name', 'vowel_name'})); 
stim_bar = bar(subject_labels, stim_acc_arr, 'FaceAlpha', .75, 'EdgeAlpha', 0);
wchance1 = yline([100/12], ':', 'Color', [0 0.4470 0.7410], 'LineWidth', chance_line_size, 'Alpha', .5 ); 
cchance1 = yline([100/4], ':', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', chance_line_size, 'Alpha', .5 ); 
vchance1 = yline([100/3], ':', 'Color', [0.9290 0.6940 0.1250],'LineWidth', chance_line_size, 'Alpha', .5 );
grid on;

stim_xtips1 = stim_bar(1).XEndPoints;
stim_ytips1 = stim_bar(1).YEndPoints;
stim_labels1 = compose("%3.1f%%", stim_bar(1).YData);
stim_text1 = text(stim_xtips1, stim_ytips1 + 3, stim_labels1,...
    'FontSize', legend_size - 3,...
    'HorizontalAlignment', 'right',...
    'VerticalAlignment', 'bottom',...
    'Rotation', 60,...
    'BackgroundColor', 'white',...
    'EdgeColor', 'black');    

stim_xtips2 = stim_bar(2).XEndPoints;
stim_ytips2 = stim_bar(2).YEndPoints;
stim_labels2 = compose("%3.1f%%", stim_bar(2).YData);
stim_text2 = text(stim_xtips2, stim_ytips2 + 3, stim_labels2,...
    'FontSize', legend_size - 3,...
    'HorizontalAlignment', 'right',...
    'VerticalAlignment', 'bottom',...
    'Rotation', 60,...
    'BackgroundColor', 'white',...
    'EdgeColor', 'black');   
    
stim_xtips3 = stim_bar(3).XEndPoints;
stim_ytips3 = stim_bar(3).YEndPoints;
stim_labels3 = compose("%3.1f%%", stim_bar(3).YData);
stim_text3 = text(stim_xtips3, stim_ytips3 + 2, stim_labels3,...
    'FontSize', legend_size - 3,...
    'HorizontalAlignment', 'left',...
    'VerticalAlignment', 'top',...
    'Rotation', 60,...
    'BackgroundColor', 'white',...
    'EdgeColor', 'black');

title('Stimulus', 'FontSize', subtitle_size);
xlabel('Subject');
ylabel('Accuracy (%)', 'FontSize', label_size);
ylim([0,80]);

inter_sub_acc_stim_ax.FontSize = label_size;
legend({plaintext(class_labels{1}), plaintext(class_labels{2}), plaintext(class_labels{3}),'Word Chance', 'Consonant Chance', 'Vowel Chance'},...
    'FontSize', legend_size);

hold off;

%%%%
inter_sub_acc_onset_ax = subplot(1, 2, 2);
onset_acc_arr = 100.*table2array(sub_onset_mean_acc_table(:,{'word_name', 'consonants_name', 'vowel_name'})); 
onset_bar = bar(subject_labels, onset_acc_arr, 'FaceAlpha', .75, 'EdgeAlpha', 0);

wchance2 = yline([100/12], ':', 'Color', [0 0.4470 0.7410], 'LineWidth', chance_line_size, 'Alpha', .5 ); 
cchance2 = yline([100/4], ':', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', chance_line_size, 'Alpha', .5 ); 
vchance2 = yline([100/3], ':', 'Color', [0.9290 0.6940 0.1250],'LineWidth', chance_line_size, 'Alpha', .5 );
grid on;

onset_xtips1 = onset_bar(1).XEndPoints;
onset_ytips1 = onset_bar(1).YEndPoints;
onset_labels1 = compose("%3.1f%%", onset_bar(1).YData);
onset_text1 = text(onset_xtips1, onset_ytips1 + 3, onset_labels1,...
    'FontSize', legend_size - 3,...
    'HorizontalAlignment', 'right',...
    'VerticalAlignment', 'bottom',...
    'Rotation', 60,...
    'BackgroundColor', 'white',...
    'EdgeColor', 'black');    

onset_xtips2 = onset_bar(2).XEndPoints;
onset_ytips2 = onset_bar(2).YEndPoints;
onset_labels2 = compose("%3.1f%%", onset_bar(2).YData);
onset_text2 = text(onset_xtips2, onset_ytips2 + 3, onset_labels2,...
    'FontSize', legend_size - 3,...
    'HorizontalAlignment', 'right',...
    'VerticalAlignment', 'bottom',...
    'Rotation', 60,...
    'BackgroundColor', 'white',...
    'EdgeColor', 'black');
    
onset_xtips3 = onset_bar(3).XEndPoints;
onset_ytips3 = onset_bar(3).YEndPoints;
onset_labels3 = compose("%3.1f%%", onset_bar(3).YData);
onset_text3 = text(onset_xtips3, onset_ytips3 + 2, onset_labels3,...
    'FontSize', legend_size - 3,...
    'HorizontalAlignment', 'left',...
    'VerticalAlignment', 'top',...
    'Rotation', 60,...
    'BackgroundColor', 'white',...
    'EdgeColor', 'black');

title('Onset', 'FontSize', subtitle_size);
xlabel('Subject');
ylabel('Accuracy (%)', 'FontSize', label_size);
ylim([0,80]);

inter_sub_acc_onset_ax.FontSize = label_size;
legend({plaintext(class_labels{1}), plaintext(class_labels{2}), plaintext(class_labels{3}),'Word Chance', 'Consonant Chance', 'Vowel Chance'},...
    'FontSize', legend_size);

figures_path = fullfile(p.grouping_path, p.current_group_value, 'Figures');
savefig(fullfile(p.figs_quickref_path, fig_quickref_filename));
exportgraphics(inter_sub_acc_fig, fullfile(p.figs_quickref_path, png_quickref_filename));
exportgraphics(inter_sub_acc_fig, fullfile(figures_path, inter_sub_acc_fig.FileName));

function s = plainformat(s)
s = regexprep(s, '[\\\^\_]','\\$0');
end
