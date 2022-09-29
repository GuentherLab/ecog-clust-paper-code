 %%% create plots comparing top coding electrodes in right vs left hemispheres
 % 
 %%% updated 2022/3/9 by AM
 

 
set(0,'DefaultFigureWindowStyle','normal')
%     set(0,'DefaultFigureWindowStyle','docked')

set_paths()

% data processed by the script load_top_encoders.m
topelc_data_filename = [ROOT_DIR, filesep, 'projectnb/busplab/Experiments/ECoG_Preprocessed_AM/topelc_data_to_surf'];



% load data to plot
load(topelc_data_filename)
nclusts = length(global_clust_num_list);
n_elcs = height(elc);

%% cluster lateralization plot
close all

plotops.ylimits = [0, 0.55];
plotops.xlimits = [0, 7];
plotops.leg_location = 'east';
plotops.yline_width = 2;
    plotops.yline_style = '--';
    plotops.yline_color =  [0.15 0.15 0.15];
plotops.annot_pos_xy =  [0.70, 0.78]; % position of pval annotation
plotops.bar_border_width =  2;
plotops.bar_colors =  [0.6 0.6 0.6; 1 1 1];
plotops.axes_line_width =  2;
plotops.axis_font_size =  13;
plotops.axes_numbers_bold =  'bold';
plotops.font =  'Arial';
plotops.fig_width_length =  [900 600];
plotops.xticklength = [0 0];
plotops.x_tick_angle =  0; % in degrees; 0=straight; positive rotates counter-clockwise
plotops.background_color =  [1 1 1]; 

hfig = figure;
hbar = bar([topelc.anytop_left_prop, 1-topelc.anytop_left_prop], 'stacked');

hfig.Color = plotops.background_color;
set(hfig,'Renderer', 'painters', 'Position', [50, 50, plotops.fig_width_length(1), plotops.fig_width_length(2)])
hbar(1).LineWidth = plotops.bar_border_width; hbar(2).LineWidth = plotops.bar_border_width; 
hbar(1).FaceColor = plotops.bar_colors(1,:);
    hbar(2).FaceColor = plotops.bar_colors(2,:);
xtickangle(gca, plotops.x_tick_angle); 
set(gca,'XTickLabels', cellstr(num2str([1:nclusts]')))
set(gca,'XColor',[0 0 0]);
set(gca,'YColor',[0 0 0]);
set(gca,'Box','off')
set(gca,'linewidth', plotops.axes_line_width)
set(gca,'FontSize', plotops.axis_font_size)
set(gca,'FontWeight', plotops. axes_numbers_bold)
set(gca,'FontName', plotops.font)
set(gca,'YTick',[0:.25:1])
xlim([plotops.xlimits]);
h=gca; h.XAxis.TickLength = plotops.xticklength;
hylabel = ylabel({'Proportion of top', 'Consonant, Vowel, or Syllable electrodes'});
    hylabel.Position = [-0.7 0.5 0.5];

 % error bars
 %%% use binomial distribution; equation from Alfonso Nieto-Castanon
 topelc.ebar_lims = NaN(nclusts,2);
 for iclust = 1:nclusts
     X = topelc.anytop_left_n(iclust); 
    N =  topelc.anytop_n(iclust); 
    alpha= 0.00001 : 0.00001 : 0.99999; 
    p=binocdf(X,N,alpha); 
    topelc.ebar_lims(iclust,:) = alpha([find(p>.975,1,'last'),find(p<.025,1,'first')]); 
 end
 hold on
 ebar_neg =  topelc.anytop_left_prop - topelc.ebar_lims(:,1); 
 ebar_pos =  -topelc.anytop_left_prop + topelc.ebar_lims(:,2); 
 h_ebar = errorbar([1:nclusts]', topelc.anytop_left_prop, ebar_neg, ebar_pos,'--');
    h_ebar.LineWidth = 0.8;
    h_ebar.LineStyle = 'none';
    h_ebar.Color = [0 0 0];
    
onset_aligned_is_left = isnan(elc.right_hemi(elc.type==1));
total_onset_proportion_left = mean(onset_aligned_is_left); % proprtion of onset aligned elcs in L hem
h_yline = yline(total_onset_proportion_left, 'LineStyle',plotops.yline_style, 'LineWidth',plotops.yline_width, 'Color',plotops.yline_color);

hleg = legend({' Left hemisphere',' Right hemisphere',''});
hleg.LineWidth = 1;
hleg.FontWeight = 'Normal';
hleg.Position =  [0.76 .89 0.23 0.07;];
hleg.EdgeColor = [1 1 1 ]; % color of border around legend

[tbl,chi2,p] = crosstab(topelc.anytop_left_n, topelc.anytop_n); 

set(gca,'LooseInset',get(gca,'TightInset')+[0.07 0 0.2 0]) % crop borders... run after positioning ylabel



%%  plot consonant vs vowel vs word, right vs left
% % % % plotops.ylimits = [0, 0.55];
plotops.xlimits = [0.2, 4.2];
plotops.leg_location = 'east';
plotops.yline_width = 2;
    plotops.yline_style = '--';
    plotops.yline_color =  [0.15 0.15 0.15];
plotops.annot_pos_xy =  [0.70, 0.78]; % position of pval annotation
plotops.bar_border_width =  2;
plotops.bar_colors =  [0.6 0.6 0.6; 1 1 1];
plotops.axes_line_width =  2;
plotops.axis_font_size =  13;
plotops.axes_numbers_bold =  'bold';
plotops.font =  'Arial';
plotops.fig_width_length =  [1000 600];
plotops.xticklength = [0 0]; 
plotops.x_tick_angle =  0; % in degrees; 0=straight; positive rotates counter-clockwise
plotops.background_color =  [1 1 1]; 

n_elc_cons_right = nnz(~isnan(elc_cons.right_hemi));
n_elc_cons_left = nnz(isnan(elc_cons.right_hemi));
n_elc_vow_right = nnz(~isnan(elc_vow.right_hemi));
n_elc_vow_left = nnz(isnan(elc_vow.right_hemi));
n_elc_word_right = nnz(~isnan(elc_word.right_hemi));
n_elc_word_left = nnz(isnan(elc_word.right_hemi));
prop_elc_cons_left = n_elc_cons_left / n_elcs; % lateralization as proportion of all electrodes
prop_elc_vow_left = n_elc_vow_left / n_elcs; % lateralization as proportion of all electrodes
prop_elc_word_left = n_elc_word_left / n_elcs; % lateralization as proportion of all electrodes

close all
hfig = figure; 
hfig.Color = plotops.background_color;
set(hfig,'Renderer', 'painters', 'Position', [50, 50, plotops.fig_width_length(1), plotops.fig_width_length(2)])

vals_to_plot = [prop_elc_cons_left, top_proportion_electrodes-prop_elc_cons_left;...
    prop_elc_vow_left, top_proportion_electrodes-prop_elc_vow_left;...
    prop_elc_word_left, top_proportion_electrodes-prop_elc_word_left];
hbar = bar(vals_to_plot, 'stacked');
set(hfig,'Renderer', 'painters', 'Position', [50, 50, plotops.fig_width_length(1), plotops.fig_width_length(2)])
hbar(1).LineWidth = plotops.bar_border_width; hbar(2).LineWidth = plotops.bar_border_width; 
hbar(1).FaceColor = plotops.bar_colors(1,:);
    hbar(2).FaceColor = plotops.bar_colors(2,:);
xtickangle(gca, plotops.x_tick_angle); 
set(gca,'XTickLabels',{'consonant', 'vowel', 'syllable'})
set(gca,'XColor',[0 0 0]);
set(gca,'YColor',[0 0 0]);
set(gca,'Box','off')
set(gca,'linewidth', plotops.axes_line_width)
set(gca,'FontSize', plotops.axis_font_size)
set(gca,'FontWeight', plotops. axes_numbers_bold)
set(gca,'FontName', plotops.font)
% set(gca,'YTick',[0:.25:1])
h=gca; h.XAxis.TickLength = plotops.xticklength;
xlim([plotops.xlimits]);
hylabel = ylabel({'Proportion of top electrodes'});
%     hylabel.Position = [-0.7 0.5 0.5];


 % error bars
 %%% use binomial distribution; equation from Alfonso Nieto-Castanon
alpha= 0.00001 : 0.00001 : 0.99999; 
n_top_elc = length(round(n_elcs*[1-top_proportion_electrodes]) : n_elcs); 
p=binocdf(n_elc_cons_left, n_top_elc, alpha); 
ebar_cons= top_proportion_electrodes * alpha([find(p>.975,1,'last'),find(p<.025,1,'first')]); 
p=binocdf(n_elc_vow_left, n_top_elc, alpha); 
ebar_vow= top_proportion_electrodes * alpha([find(p>.975,1,'last'),find(p<.025,1,'first')]); 
p=binocdf(n_elc_word_left, n_top_elc, alpha); 
ebar_word= top_proportion_electrodes * alpha([find(p>.975,1,'last'),find(p<.025,1,'first')]); 

 hold on
 ebar_neg =  [prop_elc_cons_left - ebar_cons(1); prop_elc_vow_left - ebar_vow(1); prop_elc_word_left - ebar_word(1)]; 
 ebar_pos =  [-prop_elc_cons_left + ebar_cons(2); -prop_elc_vow_left + ebar_vow(2); -prop_elc_word_left + ebar_word(2)];
 h_ebar = errorbar([1;2;3], vals_to_plot(:,1), ebar_neg, ebar_pos,'--');
    h_ebar.LineWidth = 0.8;
    h_ebar.LineStyle = 'none';
    h_ebar.Color = [0 0 0];


onset_aligned_is_left = isnan(elc.right_hemi(elc.type==1));
topelc_left_chance = mean(onset_aligned_is_left) * top_proportion_electrodes; % chance level of proportion top electrodes left
h_yline = yline(topelc_left_chance, 'LineStyle',plotops.yline_style, 'LineWidth',plotops.yline_width, 'Color',plotops.yline_color);

hleg = legend({' Left hemisphere',' Right hemisphere',''});
hleg.LineWidth = 1;
hleg.FontWeight = 'Normal';
hleg.Position =  [0.76 .89 0.25 0.07;];
hleg.EdgeColor = [1 1 1 ]; % color of border around legend

% NB: crosstab doesn't work with only 2 observations [e.g. comparing only consonant and word]
[tbl,chi2,p] = crosstab([n_elc_cons_left, n_elc_word_left]', [n_top_elc, n_top_elc]'); 





