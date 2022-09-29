function [hfig, bar_data] = plot_props(edb_nonan, stim_onset_counts, k_pct, params, per_sub_flag, plotops)

%%
field_default('plotops', 'include_unused_clusters', 0); 
field_default('plotops', 'subjects_as_colors', 0); 
field_default('plotops','add_bargraph_title',0)
    title_size = 28;
    subtitle_size = 22; 
field_default('plotops', 'ylimits',[0, 0.55])
field_default('plotops', 'xlimits',[0, 7])
field_default('plotops', 'leg_location','east')
field_default('plotops', 'yline_width', 2)
    field_default('plotops', 'yline_style', '--')
    field_default('plotops', 'yline_color', [0.15 0.15 0.15])
field_default('plotops', 'annot_pos_xy', [0.70, 0.78]); % position of pval annotation
field_default('plotops', 'bar_border_width', 2); 
field_default('plotops', 'bar_color', [0.6 0.6 0.6]);
field_default('plotops', 'axes_line_width', 2); 
field_default('plotops', 'axis_font_size', 13); 
field_default('plotops', 'axes_numbers_bold', 'bold'); 
field_default('plotops', 'font', 'Arial');
field_default('plotops', 'fig_width_length', [900 600]);
field_default('plotops', 'x_tick_angle', 0); % in degrees; 0=straight; positive rotates counter-clockwise
field_default('plotops', 'background_color', [1 1 1]); 

label_size = 24; 
legend_size = 18; 
tick_size = 18; 
stim_clust_inds = [1 2 4 5 6 8]; 
onset_clust_inds = [3:8]; 
non_stim_rows = [3, 7]; 
non_onset_rows = [1, 2];
xmark_offset = 0.1;
% align_conds = {'onset', 'stim'};
    align_conds = {'onset'}; 

close all

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

% clust_prop_fig = figure(); 
if plotops.add_bargraph_title
    if per_sub_flag
        sgtitle(sprintf('Cluster Proportions in Top %2.0f%% of %s Scores, Per Subject', k_pct * 100, class_str), 'FontSize', title_size)
    else 
        sgtitle({'Cluster Proportions in ', sprintf('Top %2.0f%% of All', k_pct * 100), sprintf('%s Scores', class_str)}, 'FontSize', title_size)
    end
end

for i_align = 1:length(align_conds)
    this_align = align_conds{i_align};
    if strcmp(this_align, 'onset')
        bar_data = gen_bar_data(onset_all_subs_counts_summary, onset_per_sub_counts, class_label);
        nanbars = onset_per_sub_counts.cluster_name(non_onset_rows);
        clustinds_this_align = onset_clust_inds;
    elseif strcmp(this_align, 'stim')
        bar_data = gen_bar_data(stim_all_subs_counts_summary, stim_per_sub_counts, class_label);
        nanbars = stim_per_sub_counts.cluster_name(non_stim_rows);
        clustinds_this_align = stim_clust_inds;
    end

    % choose which clusters to plot
    if plotops.include_unused_clusters
        clustinds_to_plot = unique([stim_clust_inds, onset_clust_inds]);
    elseif ~plotops.include_unused_clusters
        clustinds_to_plot = clustinds_this_align; 
    end
    
    % create proportions bar graph
    hfig = figure; 
    hfig.Color = plotops.background_color;
    set(hfig,'Renderer', 'painters', 'Position', [50, 50, plotops.fig_width_length(1), plotops.fig_width_length(2)])
    hold on;
    if plotops.subjects_as_colors % break out each subject's contribution to each bar
        hbar = bar(bar_data.clusters(clustinds_to_plot), bar_data.per_sub_props(clustinds_to_plot,:), 'stacked');
    elseif ~plotops.subjects_as_colors % don't differentiate subjects within a bar
        hbar = bar(sum(bar_data.per_sub_props(clustinds_to_plot,:), 2));
        set(gca,'XTickLabels',bar_data.clusters(clustinds_to_plot))
        set(gca,'XTick',1:length(clustinds_to_plot))
    end
    hbar.LineWidth = plotops.bar_border_width; 
    hbar.FaceColor = plotops.bar_color;
    xtickangle(gca, plotops.x_tick_angle); 
    set(gca,'LooseInset',get(gca,'TightInset')+[0 0 0.005 0]) % crop borders
    set(gca,'XColor',[0 0 0]);
    set(gca,'YColor',[0 0 0]);
    set(gca,'Box','off')
    set(gca,'linewidth', plotops.axes_line_width)
    set(gca,'FontSize', plotops.axis_font_size)
    set(gca,'FontWeight', plotops. axes_numbers_bold)
    set(gca,'FontName', plotops.font)
    
    % add Xs for unused clusters, if they're being plotted
    if plotops.include_unused_clusters
        scatter([nanbars(:); nanbars(:)],...
            [repmat(k_pct + xmark_offset,length(nanbars),1); 
            repmat(k_pct - xmark_offset,length(nanbars),1)],...
            'x','k', 'LineWidth', 2.5);
    end

    stim_onset_counts_top_k = bar_data.prop_data{clustinds_this_align, 2};
    counts_all = bar_data.prop_data{clustinds_this_align, 3};
    chi_pval = chitest(stim_onset_counts_top_k, counts_all, k_pct);

    h_yline = yline(k_pct, 'LineStyle',plotops.yline_style, 'LineWidth',plotops.yline_width, 'Color',plotops.yline_color);
    ylim([plotops.ylimits]);
    xlim([plotops.xlimits]);
    ylabel(sprintf('Proportion in top %2.0f%%', k_pct * 100))
    set(gca, 'FontSize', 20);
    if plotops.subjects_as_colors % break out each subject's contribution to each bar
        leg = legend(params.sub_nums_formal, 'Location', plotops.leg_location, 'FontSize', legend_size);
    end
    an = annotation('textbox', 'String', sprintf('p = %0.3f', chi_pval), 'FitBoxToText', 'on');
    an.Position(1:2) =  plotops.annot_pos_xy; % position of pval annotation
    an.FontSize = legend_size;
    
    % error bars
     %%% use binomial distribution; equation from Alfonso Nieto-Castanon
     nclusts = length(clustinds_this_align); 
     ebar_lims = NaN(nclusts,2);
     clust_top_props = sum(bar_data.per_sub_props(clustinds_to_plot,:), 2); 
     for iclust = 1:nclusts
         X = stim_onset_counts_top_k(iclust); % top elcs from this clust
        N =  counts_all(iclust); % total elcs in this clust
        alpha= 0.00001 : 0.00001 : 0.99999; 
        p=binocdf(X,N,alpha); 
        ebar_lims(iclust,:) = alpha([find(p>.975,1,'last'),find(p<.025,1,'first')]); 
     end
     hold on
     ebar_neg =  clust_top_props- ebar_lims(:,1); 
     ebar_pos =  -clust_top_props + ebar_lims(:,2); 
     h_ebar = errorbar([1:nclusts]', clust_top_props, ebar_neg, ebar_pos,'--');
        h_ebar.LineWidth = 0.8;
        h_ebar.LineStyle = 'none';
        h_ebar.Color = [0 0 0];

%     if chi_pval <= 0.05
%         an.BackgroundColor = [0.9290 0.6940 0.1250];
%     end
    hold off; 
    
end

end