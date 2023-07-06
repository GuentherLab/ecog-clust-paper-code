% generate_edb_nonan.m

close all; 

cluster_order = {'ESP-s', 'ESP-r', 'PtM-s', 'PtM-r', 'ME-sb', 'ME-sn', 'AP-r', 'AP-s'};

title_size = 28;
subtitle_size = 22; 
label_size = 24; 
legend_size = 18; 
tick_size = 18; 

edb_nonan.cluster_name = categorical(edb_nonan.cluster_name);
edb_nonan.cluster_name = reordercats(edb_nonan.cluster_name, cluster_order);

cluster_summary = groupsummary(edb_nonan, {'cluster_name', 'type'}, 'mean', {'word_accuracy_change_wo_electrode', 'cons_accuracy_change_wo_electrode', 'vowel_accuracy_change_wo_electrode'}, 'IncludeEmptyGroups', true, 'IncludeMissingGroups', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Avg Percent Change  per Cluster bar chart 
figure();
sgtitle('Average effect of removing electrodes (averaged across cluster)', 'FontSize', title_size);

subplot(2,1,1);
stim_cluster_scores = cluster_summary(cluster_summary.type == 2, :) 

bar(stim_cluster_scores{:,1}, 100 * stim_cluster_scores{:,end-2:end}); 

title('Stimulus aligned', 'FontSize', subtitle_size);
ylim([-10, 10]);
yl = ylim;
ylabel('Percent Change in Accuracy', 'FontSize', label_size);
grid on;
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',tick_size)
legend({'Word', 'Consonant', 'Vowel'}, 'Location', 'northeast', 'FontSize', legend_size)


subplot(2,1,2);
onset_cluster_scores = cluster_summary(cluster_summary.type == 1, :) 

bar(onset_cluster_scores{:,1}, 100 * onset_cluster_scores{:,end-2:end}); 

title('Voice onset aligned', 'FontSize', subtitle_size);
ylim([yl(1), yl(2)]);
ylabel('Percent Change in Accuracy', 'FontSize', label_size);
grid on;
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',tick_size)
legend({'Word', 'Consonant', 'Vowel'}, 'Location', 'northeast', 'FontSize', legend_size)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% cluster counts within accuracy distributions
figure(); 
sgtitle('Electrode counts per accuracy change', 'FontSize', title_size); 

edges = linspace(-.25, .25, 19); 

cluster_names = unique(edb_nonan.cluster_name);

plot_left = 1:2:5; 
plot_right = 2:2:6; 

for plot_row = 1:3
    class_label = p.class_names_formal{plot_row};

    vars = {'word_accuracy_change_wo_electrode',...
        'cons_accuracy_change_wo_electrode',...
        'vowel_accuracy_change_wo_electrode'};
    
    subplot(3, 2, plot_left(plot_row)); 
    histogram(discretize(edb_nonan{:, vars{plot_row}}, edges)); 
    ylabel(class_label, 'FontSize', label_size);
    if plot_row == 1
        xl = xlim; 
    else
        xlim([xl(1), xl(2)]);
    end
    
    subplot(3, 2, plot_right(plot_row));
    N_all = [];
    for cluster_idx = 1:length(cluster_names)
        cluster = cluster_names(cluster_idx);
        X = edb_nonan{edb_nonan.cluster_name == cluster, vars{plot_row}};
        [N, ~] = histcounts(X, edges);
        N_all(cluster_idx, :) = N;
    end
    b = bar(N_all', 'stacked', 'BarWidth', 1);
    legend(cluster_names, 'NumColumns', 2, 'FontSize', legend_size); 
    xlim([xl(1), xl(2)]);

end

cluster_summary = table();
for sub_num = p.sub_nums
    tt = head(edb_nonan(edb_nonan.subject == sub_num, {'subject', 'electrode', 'type', 'cluster_name', 'word_accuracy_change_wo_electrode', 'word_accuracy_change_wo_electrode'}), 4);
    cluster_summary = [cluster_summary; tt];
end

cluster_summary


function gen_fig(edb_nonan_subset)
% close all; 
title_size = 28;
subtitle_size = 22; 
label_size = 24; 
legend_size = 18; 
tick_size = 18; 

elect_list = edb_nonan.electrode; 
edges = linspace(-.25, .25, 19); 

figure();
sgtitle('Distributions of deviations from baseline accuracy', 'FontSize', title_size);
sub_num = 357; 

subplot(3,1,1);
histogram(edb_nonan.word_accuracy_change_wo_electrode(edb_nonan.subject == sub_num), edges)
title('Word', 'FontSize', subtitle_size);
xl = [-.25, .25];
yl = [0, 300];
xlim([xl(1), xl(2)]);
ylim([yl(1), yl(2)]);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',tick_size)

subplot(3,1,2);
histogram(edb_nonan.cons_accuracy_change_wo_electrode(edb_nonan.subject == sub_num), edges)
title('Consonant', 'FontSize', subtitle_size);
xlim([xl(1), xl(2)]);
ylim([yl(1), yl(2)]);
ylabel('Count', 'FontSize', label_size);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',tick_size)

subplot(3,1,3);
histogram(edb_nonan.vowel_accuracy_change_wo_electrode(edb_nonan.subject == sub_num), edges)
title('Vowel', 'FontSize', subtitle_size);
xlim([xl(1), xl(2)]);
ylim([yl(1), yl(2)]);
xlabel('\Delta Accuracy', 'FontSize', label_size);

a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',tick_size)

end
