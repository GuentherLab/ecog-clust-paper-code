 %%% create plots comparing top coding electrodes in right vs left hemispheres
 % 
 %%% updated 2022/7/30 by AM
 
 clear
 
set(0,'DefaultFigureWindowStyle','normal')
%     set(0,'DefaultFigureWindowStyle','docked')

% alpha levels for significance stars; will be adjusted for multiple comparisons
star_levels =     [inf, 0.05, 0.01, 0.001, 0.0001]; 
    star_symbols = {'', '*', '**', '***', '****'}; 

set_paths()

% data processed by the script load_top_encoders.m
topelc_data_filename = [DATA_DIR '/topelc_data_to_surf'];

% load data to plot
load(topelc_data_filename)
nclusts = length(global_clust_num_list);
n_elcs = height(elc);

%% cluster lateralization plot
close all

plotops.only_top_electrodes = 0; % if true, only plot top encoders; otherwise plot all elcs from Scott's electrode table
plotops.save_fig = 1; 
    plotops.output_resolution = 300;

plotops.inner_position = [0.17    0.1209    0.8420    0.8791];     
plotops.ylimits = [0, 1.1];
plotops.xlimits = [0, 7];
plotops.leg_location = 'east';
plotops.yline_width = 2;
    plotops.yline_style = '--';
    plotops.yline_color =  [0.15 0.15 0.15];
plotops.annot_pos_xy =  [0.70, 0.78]; % position of pval annotation
plotops.bar_border_width =  2;
plotops.bar_colors =  [0.6 0.6 0.6]; 
plotops.axes_line_width =  2;
plotops.axis_font_size =  25;
plotops.tick_font_size = 17;

plotops.axes_numbers_bold =  'bold';
plotops.font =  'Arial';
plotops.fig_width_length =  [900 600];
plotops.xticklength = [0 0];
    plotops.x_tick_angle =  0; % in degrees; 0=straight; positive rotates counter-clockwise
    plotops.xticktlabels = cellstr([repmat('C',nclusts,1), num2str([1:nclusts]')]); 
    plotops.xlabel_position = [3.6000   -0.08   -1.0000]; 
plotops.background_color =  [1 1 1]; 
plotops.ylabel_position = [-0.7 0.5 0.5];
plotops.star_font_size = 20;

plotops.nlabel_text_color = [0 0 0];
plotops.nlabel_font_size = 10; 
plotops.nlabel_font_weight = 'Bold'; % vs 'Normal'
plotops.nlabel_background_color = 'none';

if plotops.only_top_electrodes
    plotops.save_fig_filename = [fileparts(topelc_data_filename) '/figs/cluster_lateralization_bars_(top_elc)'];
    plotops.ylabel = {'Fraction of top electrodes','in left hemisphere'};
    n_per_clust = topelc.anytop_n; 
    left_n_per_clust = topelc.anytop_left_n; 
    % for topelc only, set yline =  L proportion of topelc only
    proportion_left_overall = sum(topelc.anytop_left_n) / sum(topelc.anytop_n);
    labels_n_elc = [topelc.anytop_left_n, topelc.anytop_n]; 
elseif ~plotops.only_top_electrodes
    plotops.save_fig_filename = [fileparts(topelc_data_filename) '/figs/cluster_lateralization_bars_(all_elc)'];
    plotops.ylabel = {'Fraction of responsive electrodes','in left hemisphere'};
    n_per_clust = topelc.n_elc; 
    left_n_per_clust = topelc.n_elc_RL(:,2); 
    % set yline at L proportion of all electrodes
    onset_aligned_is_left = isnan(elc.right_hemi(elc.type==1));
    proportion_left_overall = mean(onset_aligned_is_left); % proprtion of onset aligned elcs in L hem
    labels_n_elc = [topelc.n_elc_RL(:,2), topelc.n_elc]; 
end
left_prop_by_clust = left_n_per_clust ./ n_per_clust; 

hfig = figure;
hbar = bar([left_prop_by_clust], 'stacked');

hfig.Color = plotops.background_color;
set(hfig,'Renderer', 'painters', 'Position', [50, 50, plotops.fig_width_length(1), plotops.fig_width_length(2)])
hbar(1).LineWidth = plotops.bar_border_width;
hbar(1).FaceColor = plotops.bar_colors(1,:);

hax = gca;
hax.InnerPosition = plotops.inner_position;
xtickangle(hax, plotops.x_tick_angle); 
hax.FontSize = plotops.tick_font_size;
set(gca,'XTickLabels', plotops.xticktlabels)
set(gca,'XColor',[0 0 0]);
set(gca,'YColor',[0 0 0]);
set(gca,'Box','off')
set(gca,'linewidth', plotops.axes_line_width)
set(gca,'FontWeight', plotops. axes_numbers_bold)
set(gca,'FontName', plotops.font)
set(gca,'YTick',[0:.25:1])
xlim(plotops.xlimits);
ylim(plotops.ylimits); 
hax.XAxis.TickLength = plotops.xticklength;
hxlabel = xlabel('Cluster');
    hxlabel.Position = plotops.xlabel_position;
    hxlabel.FontSize = plotops.axis_font_size;
hylabel = ylabel(plotops.ylabel);
    hylabel.Position = plotops.ylabel_position;
    hylabel.FontSize = plotops.axis_font_size;

% adjust significance stars for Bonferroni correction
adjusted_star_levels = star_levels / nclusts;
alpha_bonf = 0.05 / nclusts; % 2-tailed alpha = 0.05 w/ Bonferroni correction

 % error bars
 %%% use binomial distribution; equation from Alfonso Nieto-Castanon
ebar_lims = NaN(nclusts,2);
p_binotest = NaN(nclusts,1); 
q_binotest = NaN(nclusts,1); 
sgn_stars = cell(nclusts,1);
sgn_stars_bonf = cell(nclusts,1);
 for iclust = 1:nclusts
     X = left_n_per_clust(iclust); 
    N =  n_per_clust(iclust); 
    alpha= 0.00001 : 0.00001 : 0.99999; 
    p=binocdf(X,N,alpha); 
    ebar_lims(iclust,:) = alpha([find(p>.975,1,'last'),find(p<.025,1,'first')]); 
    % use X-1 because matlab binocdf with ‘upper’ computes the probability of getting above X...
    %       .... whereas we want to ask (in the upper half of the test) the inclusive probably of getting our outcome or above
    p_binotest(iclust) = 2 * min([binocdf(X,N,proportion_left_overall), binocdf(X-1,N,proportion_left_overall,'upper')]); % 2* because 2-tailed
    starind_bonf = find(p_binotest(iclust) < adjusted_star_levels, 1, 'last');
    sgn_stars_bonf{iclust} = star_symbols{starind_bonf}; 
 end
 sgn_lat_bonf = p_binotest < alpha_bonf; % significance after bonf correction
 
 % FDR correction
 q_binotest = conn_fdr(p_binotest); 
 for iclust = 1:nclusts % need to run FDR and set star values after computing all base p-values
     starind = find(q_binotest(iclust) < star_levels, 1, 'last');
     sgn_stars{iclust} = star_symbols{starind};
 end
 sgn_lat_fdr = q_binotest < 0.05;
 
lateralization =  table(left_prop_by_clust, n_per_clust, left_n_per_clust, p_binotest, q_binotest, sgn_lat_fdr, sgn_stars, sgn_lat_bonf, sgn_stars_bonf)
 
hold on
 ebar_neg =  left_prop_by_clust - ebar_lims(:,1); 
 ebar_pos =  -left_prop_by_clust + ebar_lims(:,2); 
 h_ebar = errorbar([1:nclusts]', left_prop_by_clust, ebar_neg, ebar_pos,'--');
    h_ebar.LineWidth = 0.8;
    h_ebar.LineStyle = 'none';
    h_ebar.Color = [0 0 0];
    
h_yline = yline(proportion_left_overall, 'LineStyle',plotops.yline_style, 'LineWidth',plotops.yline_width, 'Color',plotops.yline_color);

star_scat = textscatter(1:nclusts, 1.05*ones(1,nclusts), sgn_stars, 'FontSize',plotops.star_font_size);

% % % % % % % % % hleg = legend({' Left hemisphere',' Right hemisphere',''});
% % % % % % % % % hleg.LineWidth = 1;
% % % % % % % % % hleg.FontWeight = 'Normal';
% % % % % % % % % hleg.Position =  [0.76 .89 0.23 0.07;];
% % % % % % % % % hleg.EdgeColor = [1 1 1 ]; % color of border around legend

% chi-test against the null hypothesis that left hem electrodes are evenly distributed acrocss clusters
[chi_significant, chi_p, chi_stats] = chi2gof([1:nclusts]', 'Frequency',left_n_per_clust, 'Expected',proportion_left_overall*n_per_clust, 'Emin',0)

set(gca,'LooseInset',get(gca,'TightInset')+[0.1 0 0.0 0]) % crop borders... run after positioning ylabel

nlabel_text = strrep(cellstr([num2str(labels_n_elc(:,1)), repmat('/',nclusts,1), strtrim(num2str(labels_n_elc(:,2)))]),' ','');
nlabel_yval = min(ebar_lims(:)) / 2; % labels halfway up to bottom of lowest errorbar
nlabel_scat = textscatter(1:nclusts, nlabel_yval*ones(1,nclusts), nlabel_text,...
    'FontSize',plotops.nlabel_font_size, 'ColorData',plotops.nlabel_text_color,...
    'FontWeight',plotops.nlabel_font_weight, 'BackgroundColor',plotops.nlabel_background_color);


% save fig
if  plotops.save_fig
    print(hfig, plotops.save_fig_filename, '-dpng', ['-r', num2str(plotops.output_resolution)]); % save the image to file
end

%%  plot consonant vs vowel vs word, right vs left
clear plotops

 plotops.save_fig = 1; 
    plotops.save_fig_filename = [fileparts(topelc_data_filename) '/figs/cons_vow_syl_lateralization_bars'];
    plotops.output_resolution = 300;

plotops.ylimits = [0, 1];
plotops.xlimits = [0.2, 3.7];
plotops.leg_location = 'east';
plotops.yline_width = 2;
    plotops.yline_style = '--';
    plotops.yline_color =  [0.15 0.15 0.15];
plotops.annot_pos_xy =  [0.70, 0.78]; % position of pval annotation
plotops.bar_border_width =  2;

% plotops.bar_colors =  [1 0 0; ...    % RGB bars
%                       0 1 0;...
%                       0 0 1]; 
plotops.bar_colors =  [0.6 0.6 0.6; ...  % gray bars
                      0.6 0.6 0.6;...
                      0.6 0.6 0.6];                   

plotops.axes_line_width =  2;
plotops.axis_font_size =  13;
plotops.axes_numbers_bold =  'bold';
plotops.font =  'Arial';
plotops.fig_width_length =  [900 600];
plotops.xticklength = [0 0]; 
plotops.x_tick_angle =  0; % in degrees; 0=straight; positive rotates counter-clockwise
plotops.background_color =  [1 1 1]; 
plotops.star_font_size = 20;

plotops.nlabel_text_color = [0 0 0];
plotops.nlabel_font_size = 12; 
plotops.nlabel_font_weight = 'Bold'; % vs 'Normal'
plotops.nlabel_background_color = [1 1 1];

feats = {'cons', 'vow', 'word'};
    plotops.xticklabels = {'Consonant', 'Vowel', 'Syllable'};
    nfeats = length(feats);
    
n_elc_cons_right = nnz(~isnan(elc_cons.right_hemi));
n_elc_cons_left = nnz(isnan(elc_cons.right_hemi));
n_elc_vow_right = nnz(~isnan(elc_vow.right_hemi));
n_elc_vow_left = nnz(isnan(elc_vow.right_hemi));
n_elc_syl_right = nnz(~isnan(elc_word.right_hemi));
n_elc_syl_left = nnz(isnan(elc_word.right_hemi));

nans = nan(nfeats,1);
cel = cell(nfeats,1);
featlat = table(feats',plotops.xticklabels',nans,nans,nans,nans,nans,cel,'VariableNames',...
        {'feat','yticklabel','left_prop','n_total','n_left','p_binotest','sgn_lateralization','sgn_stars'});
featlat.q_binotest = nans;
featlat.sgn_stars_bonf = cel; 
featlat.n_total = [n_elc_cons_right + n_elc_cons_left; n_elc_vow_right + n_elc_vow_left; n_elc_syl_right + n_elc_syl_left]; 
featlat.n_left = [n_elc_cons_left; n_elc_vow_left; n_elc_syl_left]; 
featlat.left_prop = featlat.n_left ./ featlat.n_total;

% set yline at L proportion of all electrodes
onset_aligned_is_left = isnan(elc.right_hemi(elc.type==1));
proportion_left_overall = mean(onset_aligned_is_left); % proprtion of onset aligned elcs in L hem

% adjust significance stars for 2-tailed test, then for Bonferroni correction
adjusted_star_levels = star_levels / nfeats;
alpha_bonf = 0.05 / nfeats; % alpha = 0.05 w/ Bonferroni correction

close all
hfig = figure; 
hfig.Color = plotops.background_color;
set(hfig,'Renderer', 'painters', 'Position', [50, 50, plotops.fig_width_length(1), plotops.fig_width_length(2)])

hbar = bar(featlat.left_prop);
hbar.LineWidth = plotops.bar_border_width;
hbar.FaceColor = 'flat';

 for ifeat = 1:nfeats
    thisfeat = feats{ifeat};
    X = featlat.n_left(ifeat); % n left elcs which are top coders of this type
    N = featlat.n_total(ifeat); % total top coders of this type
    alpha= 0.00001 : 0.00001 : 0.99999; 
    p=binocdf(X,N,alpha); 
    featlat.ebar_lims(ifeat,:) = alpha([find(p>.975,1,'last'),find(p<.025,1,'first')]); 
    % use X-1 because matlab binocdf with ‘upper’ computes the probability of getting above X...
    %       .... whereas we want to ask (in the upper half of the test) the inclusive probably of getting our outcome or above
    featlat.p_binotest(ifeat) = 2 * min([binocdf(X,N,proportion_left_overall), binocdf(X-1,N,proportion_left_overall,'upper')]);  % 2* because 2-tailed
    starind_bonf = find(featlat.p_binotest(ifeat) < adjusted_star_levels, 1, 'last');
    featlat.sgn_stars_bonf{ifeat} = star_symbols{starind_bonf}; 
    
    % bar colors
    hbar.CData(ifeat,:) = plotops.bar_colors(ifeat,:);
 end
featlat.sgn_lat_bonf = featlat.p_binotest < alpha_bonf;

 % FDR correction
 featlat.q_binotest = conn_fdr(featlat.p_binotest); 
 for ifeat = 1:nfeats % need to run FDR and set star values after computing all base p-values
     starind = find(featlat.q_binotest(ifeat) < star_levels, 1, 'last');
     featlat.sgn_stars{ifeat} = star_symbols{starind};
 end


set(hfig,'Renderer', 'painters', 'Position', [50, 50, plotops.fig_width_length(1), plotops.fig_width_length(2)])
xtickangle(gca, plotops.x_tick_angle); 
hax = gca; 
set(gca,'XTickLabels',plotops.xticklabels)
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
ylim(plotops.ylimits); 
hylabel = ylabel({'Fraction of top electrodes','in left hemisphere'});
    hylabel.Position = [-0.05   .5   0];

 hold on
 ebar_neg =  featlat.left_prop - featlat.ebar_lims(:,1); 
 ebar_pos =  -featlat.left_prop + featlat.ebar_lims(:,2); 
 h_ebar = errorbar([1:nfeats]', featlat.left_prop, ebar_neg, ebar_pos,'--');
    h_ebar.LineWidth = 0.8;
    h_ebar.LineStyle = 'none';
    h_ebar.Color = [0 0 0];

h_yline = yline(proportion_left_overall, 'LineStyle',plotops.yline_style, 'LineWidth',plotops.yline_width, 'Color',plotops.yline_color);

%%%% significance stars
star_ybuffer = 0.1; 
star_yval = star_ybuffer + max(featlat.ebar_lims(:)); 
star_scat = textscatter(1:nfeats, star_yval*ones(1,nfeats), featlat.sgn_stars, 'FontSize',plotops.star_font_size);

% labels for n electrodes per group
labels_n_elc = [featlat.n_left, featlat.n_total]; 
nlabel_text = strrep(cellstr([num2str(labels_n_elc(:,1)), repmat('/',nfeats,1), strtrim(num2str(labels_n_elc(:,2)))]),' ','');
nlabel_yval = min(featlat.ebar_lims(:)) / 2; % labels halfway up to bottom of lowest errorbar
nlabel_scat = textscatter(1:nfeats, nlabel_yval*ones(1,nfeats), nlabel_text,...
    'FontSize',plotops.nlabel_font_size, 'ColorData',plotops.nlabel_text_color,...
    'FontWeight',plotops.nlabel_font_weight, 'BackgroundColor',plotops.nlabel_background_color);

% chi-test against the null hypothesis that left hem electrodes are evenly distributed acrocss clusters
[chi_significant, chi_p, chi_stats] = chi2gof([1:nfeats]', 'Frequency',featlat.n_left, 'Expected',proportion_left_overall*featlat.n_total, 'Emin',0)

% % % % % % % % % % hleg = legend({' Left hemisphere',' Right hemisphere',''});
% % % % % % % % % % hleg.LineWidth = 1;
% % % % % % % % % % hleg.FontWeight = 'Normal';
% % % % % % % % % % hleg.Position =  [0.76 .89 0.25 0.07;];
% % % % % % % % % % hleg.EdgeColor = [1 1 1 ]; % color of border around legend


% save fig
if  plotops.save_fig
    print(hfig, plotops.save_fig_filename, '-dpng', ['-r', num2str(plotops.output_resolution)]); % save the image to file
end



