 %%% create plots comparing top coding electrodes in right vs left hemispheres
 % 
 %%% updated 2022/3/9 by AM
 

 
 reload_electrode_data = 0; % takes a few minutes... turn off if data already loaded
 load_preloaded_electrode_data = 1; % alternative to running the slow loading function... ignored if reload_electrode_data == 1
    preloaded_electrode_data_filename = '/usr2/postdoc/amsmeier/leave_one_out_data';
 
 addpath('/usr2/postdoc/amsmeier')
if  reload_electrode_data
    
     load_top_coding_electrodes 
elseif load_preloaded_electrode_data
    load(preloaded_electrode_data_filename,'ow_sorted') % onset data sorted by word coding
end
if exist('ow_sorted','var') % only analyze onset data
    elc = sortrows(ow_sorted,{'subject','electrode'}); % sort electrodes by subject and electrode ID
end
elc.cluster_name = string(elc.cluster_name);
clearvars -except elc

%%
 top_proportion_electrodes = 0.3; % use this proportion of the electrodes as 'top coders'

 %%%%%%%%%%%%%%%% organize data %%%%%%%%%%
n_elcs = height(elc);  
falsecol = false(n_elcs,1);
elc = [elc, table(falsecol, falsecol, falsecol, falsecol, 'VariableNames', {'top_cons_coder','top_vow_coder','top_word_coder','top_cons_or_word'})];                              

% create table with top-30% coding electrodes data by cluster
clustlist = {'PtM-s','PtM-r','ME-sb','ME-sn','AP-r','AP-s'}';
nclusts = length(clustlist); 

% make tables of top coders
top_inds = round(n_elcs*[1-top_proportion_electrodes]) : n_elcs; 
n_top_elc = length(top_inds);
elc = sortrows(elc,'cons_accuracy_change_wo_electrode');
    elc.top_cons_coder(top_inds) = true; 
    elc_cons = elc(top_inds,:); 
elc = sortrows(elc,'vowel_accuracy_change_wo_electrode');
    elc.top_vow_coder(top_inds) = true; 
    elc_vow = elc(top_inds,:); 
elc = sortrows(elc,'word_accuracy_change_wo_electrode');
    elc.top_word_coder(top_inds) = true; 
    elc_word = elc(top_inds,:); 
    
elc.top_cons_or_word(elc.top_cons_coder | elc.top_word_coder) = true; 
    elc_cons_or_word = elc(elc.top_cons_or_word,:); 
    elc = sortrows(elc,{'subject','electrode'}); % sort electrodes by subject and electrode ID
    
nancol = NaN(nclusts, 1);
nan2 = [nancol, nancol]; 
topelc = table(clustlist, nancol,   nancol,              nan2,       nancol,              nan2,         nancol,            nan2,       nancol,          nancol,  'VariableNames',...
                    {'clust',       'n_elc', 'cons_left_prop', 'cons_RL', 'vow_left_prop', 'vow_RL', 'word_left_prop', 'word_RL', 'cons_word_left_n', 'cons_word_left_prop'}); 
    
for iclust = 1:nclusts
    thisclust = clustlist{iclust};
    clustrows = strcmp(elc.cluster_name, thisclust);
    topelc.n_elc(iclust) = nnz(clustrows);
    
    clustrows_top_cons = clustrows & elc.top_cons_coder; 
    n_top_cons_this_clust = nnz(clustrows_top_cons); % number of electrodes in this cluster which are top cons coderclc
    
    n_top_cons_this_clust_leftside = nnz(clustrows_top_cons & isnan(elc.right_hemi)); % number electrodes in this clust which are top coders and in L hem 
    topelc.cons_RL(iclust,:) = [n_top_cons_this_clust-n_top_cons_this_clust_leftside  , n_top_cons_this_clust_leftside]; % number of R and L electrodes
    topelc.cons_left_prop(iclust) = n_top_cons_this_clust_leftside / n_top_cons_this_clust; 
    
    clustrows_top_vow = clustrows & elc.top_vow_coder; 
    n_top_vow_this_clust = nnz(clustrows_top_vow); % number of electrodes in this cluster which are top coder
    n_top_vow_this_clust_leftside = nnz(clustrows_top_vow & isnan(elc.right_hemi)); % number electrodes in this clust which are top coders and in L hem 
    topelc.vow_RL(iclust,:) = [n_top_vow_this_clust-n_top_vow_this_clust_leftside  , n_top_vow_this_clust_leftside]; % number of R and L electrodes
    topelc.vow_left_prop(iclust) = n_top_vow_this_clust_leftside / n_top_vow_this_clust; 
    
    clustrows_top_word = clustrows & elc.top_word_coder; 
    n_top_word_this_clust = nnz(clustrows_top_word); % number of electrodes in this cluster which are top coder
    n_top_word_this_clust_leftside = nnz(clustrows_top_word & isnan(elc.right_hemi));  
    topelc.word_RL(iclust,:) = [n_top_word_this_clust-n_top_word_this_clust_leftside  , n_top_word_this_clust_leftside]; % number of R and L electrodes
    topelc.word_left_prop(iclust) = n_top_word_this_clust_leftside / n_top_word_this_clust; 
    
    clustrows_top_cons_or_word = clustrows & elc.top_cons_or_word; 
    n_top_cons_or_word_this_clust = nnz(clustrows_top_cons_or_word); % number of electrodes in this cluster which are top coder
    topelc.cons_word_n(iclust,1) = n_top_cons_or_word_this_clust; 
    n_top_cons_or_word_this_clust_leftside = nnz(clustrows_top_cons_or_word & isnan(elc.right_hemi));
    topelc.cons_word_left_n(iclust) = n_top_cons_or_word_this_clust_leftside; 
    topelc.cons_word_left_prop(iclust) = n_top_cons_or_word_this_clust_leftside / n_top_cons_or_word_this_clust; 
end

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
hbar = bar([topelc.cons_word_left_prop, 1-topelc.cons_word_left_prop], 'stacked');

hfig.Color = plotops.background_color;
set(hfig,'Renderer', 'painters', 'Position', [50, 50, plotops.fig_width_length(1), plotops.fig_width_length(2)])
hbar(1).LineWidth = plotops.bar_border_width; hbar(2).LineWidth = plotops.bar_border_width; 
hbar(1).FaceColor = plotops.bar_colors(1,:);
    hbar(2).FaceColor = plotops.bar_colors(2,:);
xtickangle(gca, plotops.x_tick_angle); 
set(gca,'XTickLabels', clustlist)
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
hylabel = ylabel({'Proportion of top', 'Word or Consonant electrodes'});
    hylabel.Position = [-0.7 0.5 0.5];

 % error bars
 %%% use binomial distribution; equation from Alfonso Nieto-Castanon
 topelc.ebar_lims = NaN(nclusts,2);
 for iclust = 1:nclusts
     X = topelc.cons_word_left_n(iclust); 
    N =  topelc.cons_word_n(iclust); 
    alpha= 0.00001 : 0.00001 : 0.99999; 
    p=binocdf(X,N,alpha); 
    topelc.ebar_lims(iclust,:) = alpha([find(p>.975,1,'last'),find(p<.025,1,'first')]); 
 end
 hold on
 ebar_neg =  topelc.cons_word_left_prop - topelc.ebar_lims(:,1); 
 ebar_pos =  -topelc.cons_word_left_prop + topelc.ebar_lims(:,2); 
 h_ebar = errorbar([1:nclusts]', topelc.cons_word_left_prop, ebar_neg, ebar_pos,'--');
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

[tbl,chi2,p] = crosstab(topelc.cons_word_left_n, topelc.cons_word_n); 

set(gca,'LooseInset',get(gca,'TightInset')+[0.07 0 0.2 0]) % crop borders... run after positioning ylabel



%%  plot consonant vs vowel, right vs left
% % % % plotops.ylimits = [0, 0.55];
plotops.xlimits = [0.2, 3];
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
n_elc_word_right = nnz(~isnan(elc_word.right_hemi));
n_elc_word_left = nnz(isnan(elc_word.right_hemi));
prop_elc_cons_left = n_elc_cons_left / n_elcs; % lateralization as proportion of all electrodes
prop_elc_word_left = n_elc_word_left / n_elcs; % lateralization as proportion of all electrodes

close all
hfig = figure; 
hfig.Color = plotops.background_color;
set(hfig,'Renderer', 'painters', 'Position', [50, 50, plotops.fig_width_length(1), plotops.fig_width_length(2)])

vals_to_plot = [prop_elc_cons_left, top_proportion_electrodes-prop_elc_cons_left;...
    prop_elc_word_left, top_proportion_electrodes-prop_elc_word_left];
hbar = bar(vals_to_plot, 'stacked');
set(hfig,'Renderer', 'painters', 'Position', [50, 50, plotops.fig_width_length(1), plotops.fig_width_length(2)])
hbar(1).LineWidth = plotops.bar_border_width; hbar(2).LineWidth = plotops.bar_border_width; 
hbar(1).FaceColor = plotops.bar_colors(1,:);
    hbar(2).FaceColor = plotops.bar_colors(2,:);
xtickangle(gca, plotops.x_tick_angle); 
set(gca,'XTickLabels',{'consonant', 'word'})
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
hylabel = ylabel({'Proportion of top', 'Word or Consonant electrodes'});
%     hylabel.Position = [-0.7 0.5 0.5];


 % error bars
 %%% use binomial distribution; equation from Alfonso Nieto-Castanon
alpha= 0.00001 : 0.00001 : 0.99999; 
p=binocdf(n_elc_cons_left, n_top_elc, alpha); 
ebar_cons= top_proportion_electrodes * alpha([find(p>.975,1,'last'),find(p<.025,1,'first')]); 
p=binocdf(n_elc_word_left, n_top_elc, alpha); 
ebar_word= top_proportion_electrodes * alpha([find(p>.975,1,'last'),find(p<.025,1,'first')]); 

 hold on
 ebar_neg =  [prop_elc_cons_left - ebar_cons(1); prop_elc_word_left - ebar_word(1)]; 
 ebar_pos =  [-prop_elc_cons_left + ebar_cons(2); -prop_elc_word_left + ebar_word(2)];
 h_ebar = errorbar([1;2], vals_to_plot(1,:), ebar_neg, ebar_pos,'--');
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

% crosstab doesn't work with only 2 observations
% % [tbl,chi2,p] = crosstab([n_elc_cons_left, n_elc_word_left]', [n_top_elc, n_top_elc]'); 





