 %%%% scatterplot top encoding proportion against a continuous cluster variable
 %.... including errorbars
 %.... then perform logistic regression
clear
 close all
 
set_paths()
topelc_data_filename = [DATA_DIR '/topelc_data_to_surf'];
load(topelc_data_filename) % load data to plot

%% parameters
xscalefact = 0.001; % convert predictor units by multiplying by this factor

plotops.predictor_var = 'width';
%     plotops.precictor_var = 'width_to_peak';
%     plotops.precictor_var = 'width_from_peak';
% plotops.precictor_var = 'start';
%     plotops.precictor_var = 'peak';
%     plotops.precictor_var = 'end';
% plotops.precictor_var = 'onset';
%     plotops.precictor_var = 'offset';


plotops.save_fig_filename = [fileparts(topelc_data_filename) '/figs/encoding_vs_', plotops.predictor_var];
    plotops.save_pdf = 0; % vector
    plotops.save_tif = 1; % raster
        plotops.output_resolution = 300; % only for raster
plotops.ylimits = [0, 1];
plotops.xlimits = xscalefact * [500, 1950];
plotops.leg_location = 'east';
plotops.background_color =  [1 1 1]; 
plotops.fig_x_y_width_height =  [50 50 1600 600]; % y coord may have to be changed per computer

% plotops.annot_pos_xy =  [0.70, 0.78]; % position of pval annotation
plotops.axes_line_width =  2;
plotops.axis_font_size =  20;
plotops.tick_font_size = 13; 
plotops.axes_numbers_bold =  'bold';
plotops.font =  'Arial';
plotops.box_on_off = 'off';
plotops.background_color =  [1 1 1]; 
plotops.padding = 'normal'; 

plotops.ylabel_position_xy = [xscalefact*260, 0.45];
plotops.yticks = [0:.25:1]; 

plotops.xlabel_position_xy  = [xscalefact*1225,   -0.07];
plotops.xticklength = [0 0];
plotops.x_tick_angle =  0; % in degrees; 0=straight; positive rotates counter-clockwise
% plotops.xticktlabels = num2str(1:6);
plotops.x_ax_label = 'Response width (sec)';

% scatterplot parameters
plotops.text_to_scat = cellstr([repmat('C',nclusts,1), num2str([1:nclusts]')]); % clust names for article
%   plotops.text_to_scat = topelc.clust; %%% old cluster names
xoffset = -60;
    plotops.text_scat_x_offset = xoffset*ones(6,3); % shift text labels rightward this much relative to markers
%     plotops.text_scat_x_offset(1,3) = -xoffset; % flip to other side
    plotops.text_scat_x_offset(2,2) = -xoffset; % flip to other side
    plotops.text_scat_x_offset(3,:) = -xoffset; % flip to other side
    plotops.text_scat_x_offset(4,2:3) = -xoffset; % flip to other side|
%     plotops.text_scat_x_offset(5,:) = -xoffset; % flip to other side
yoffset = -0.025;
    plotops.text_scat_y_offset = yoffset *ones(6,3); % shift text labels up this much relative to markers
    plotops.text_scat_y_offset(1,3) = -yoffset; % flip to other side
    plotops.text_scat_y_offset(2,1:2) = -yoffset; % flip to other side
    plotops.text_scat_y_offset(3,2:3) = -yoffset; % flip to other side
    plotops.text_scat_y_offset(4,2:3) = -yoffset; % flip to other side
    plotops.text_scat_y_offset(5,1:2) = -yoffset; % flip to other side
plotops.text_fontweight = 'bold';
plotops.text_scat_fontsize = 10; 
plotops.marker_size = 42; 
plotops.marker_color = [0 0 0];
plotops.text_background_color = [1 1 1];    % RBG or 'none'

plotops.title_fontsize = 17; 
plotops.title_fontweight = 'bold';
%     plotops.title_fontweight = 'normal';    

% plotops.yline_width = 2;
%     plotops.yline_style = '--';
%     plotops.yline_color =  [0.15 0.15 0.15];

plotops.trend_width = 2; 
plotops.trend_color = [0.4 0.4 1]; 
plotops.trend_style = '-';

plotops.ebar_line_width = 0.5;
plotops.ebar_line_style = 'none'; % connecting lines
plotops.ebar_color = [0 0 0];
plotops.ebar_vline_style = 'dotted';  %% vertical line style; options: 'solid' | 'dashed' | 'dotted' | 'dashdot' | 'none'

feats = {'cons', 'vow', 'word'};
    featlabels = {'Consonant', 'Vowel', 'Syllable'};
    nfeats = length(feats);

    set(0,'DefaultFigureWindowStyle','normal')
%     set(0,'DefaultFigureWindowStyle','docked')
    
    
%%% organize data
n_elcs = height(elc);

% assign predictor var values to electrode table
for iclust = 1:nclusts
    matchrows = strcmp(elc.cluster_name, topelc.clust{iclust});
    elc{matchrows, plotops.predictor_var} = topelc{iclust, plotops.predictor_var};
end

%%% scatterplots
close all
hfig = figure;
hfig.Color = plotops.background_color;
set(hfig,'Renderer', 'painters', 'Position', [plotops.fig_x_y_width_height ])
htile = tiledlayout(1,nfeats, 'Padding',plotops.padding);
stats = table;
    stats.feat = featlabels'; 
for ifeat = 1:nfeats
    thisfeat = feats{ifeat}; 
        thisfeatlabel = featlabels{ifeat}; 
    nexttile

    hold on
    hscat(ifeat) = scatter(xscalefact * topelc{:,plotops.predictor_var}, topelc{:,[thisfeat,'_prop']});
        hscat(ifeat).MarkerFaceColor = plotops.marker_color;
        hscat(ifeat).MarkerEdgeColor = plotops.marker_color;
        hscat(ifeat).SizeData = plotops.marker_size; 
        
    hax(ifeat) = gca; 
    hax(ifeat).XTickLabelRotation = plotops.x_tick_angle; 
    set(gca,'Box',plotops.box_on_off)
    set(gca,'linewidth', plotops.axes_line_width)
    hax(ifeat).FontSize = plotops.tick_font_size;
    set(gca,'FontWeight', plotops. axes_numbers_bold)
    set(gca,'FontName', plotops.font)
    set(gca,'YTick',plotops.yticks)
    xlim([plotops.xlimits]);
    ylim(plotops.ylimits); 
    h=gca; h.XAxis.TickLength = plotops.xticklength;
    yticklabs = hax(ifeat).YTick; % save for later
    hax(ifeat).YTickLabels = '';

%     xlab = [upper(plotops.predictor_var(1)), plotops.predictor_var(2:end)]; % capitalize
    xlab = plotops.x_ax_label; 
    hxlabel(ifeat) = xlabel(xlab); 
        hxlabel(ifeat).Position = plotops.xlabel_position_xy;
        hxlabel(ifeat).FontSize = plotops.axis_font_size;
%     hylabel(ifeat) = ylabel({['Proportion of top ',thisfeatlabel,'-encoding electrodes']});
%         hylabel(ifeat).Position = plotops.ylabel_position;


   % logistic regression
   mdl = fitglm(xscalefact * elc{:,plotops.predictor_var}, categorical(elc{:,['top_',thisfeat,'_coder']}), 'Distribution','binomial');
   stats.pval(ifeat,1:2) = mdl.Coefficients.pValue';
   stats.beta(ifeat,1:2) = [mdl.Coefficients{'(Intercept)','Estimate'}, mdl.Coefficients{'x1','Estimate'}];
   stats.odds_rat(ifeat) = exp(stats.beta(ifeat,2)); 
   stats.loglike(ifeat) = mdl.LogLikelihood;

    %%% trendline
            %linear fit, weighted by number of electrodes in each cluster
            %     [trend_p,trend_S,trend_mu] = polyfit(xscalefact * elc{:,plotops.predictor_var},elc{:,['top_',thisfeat,'_coder']},1);
            %     trend_x = xlim;
            %     trend_y = polyval(trend_p,trend_x,trend_S,trend_mu);
    % logreg fit
    trend_x = linspace(plotops.xlimits(1), plotops.xlimits(2), 1000)';
    trend_y = predict(mdl, trend_x);
    htrend = plot(trend_x,trend_y);
        htrend.LineWidth = plotops.trend_width; 
        htrend.Color = plotops.trend_color; 
        htrend.LineStyle = plotops.trend_style;

    %%% error bars
     ebar_neg =  topelc{:,[thisfeat,'_prop']} - topelc{:,[thisfeat,'_ebar']}(:,1); 
     ebar_pos =  -topelc{:,[thisfeat,'_prop']} + topelc{:,[thisfeat,'_ebar']}(:,2); 
     h_ebar(ifeat) = errorbar(xscalefact * topelc{:,plotops.predictor_var}, topelc{:,[thisfeat,'_prop']}, ebar_neg, ebar_pos,'--');
        h_ebar(ifeat).LineWidth = plotops.ebar_line_width;
        h_ebar(ifeat).LineStyle = plotops.ebar_line_style;
        h_ebar(ifeat).Color = plotops.ebar_color;
        h_ebar(ifeat).Bar.LineStyle = plotops.ebar_vline_style;
        
    htitle(ifeat) = title(thisfeatlabel);
        htitle(ifeat).FontSize = plotops.title_fontsize; 
        htitle(ifeat).FontWeight = plotops.title_fontweight;
    
   % text scatter last to put text background on top of other items
   tscat(ifeat) = textscatter(xscalefact * [topelc{:,plotops.predictor_var} + plotops.text_scat_x_offset(:,ifeat)], ...
            topelc{:,[thisfeat,'_prop']} + plotops.text_scat_y_offset(:,ifeat), plotops.text_to_scat);
        tscat(ifeat).TextDensityPercentage = 100; % don't replace text with dots
        tscat(ifeat).FontName = plotops.font;
        tscat(ifeat).FontSize = plotops.text_scat_fontsize;
        tscat(ifeat).FontWeight = plotops.text_fontweight; 
        tscat(ifeat).BackgroundColor = plotops.text_background_color; 
        

end
hax(1).YTickLabels = yticklabs;
hax(1).YLabel.String = 'Fraction of top electrodes';
    hax(1).YLabel.Position(1:2) = plotops.ylabel_position_xy;
    hax(1).YLabel.FontSize = plotops.axis_font_size;


stats.qval = conn_fdr(stats.pval); % FDR correction
stats
    
    
    
% stats = cellfun(@struct,mdl'); 
% stats = struct2table(rmfield(stats,{'Variables','VariableNames'})); 
    
% stats = struct2table(stats); 
% stats.beta = [cell2mat(stats.beta')]';
% stats.p = [cell2mat(stats.p')]';
% 
% stats.feat = feats';
% stats = movevars(stats,{'feat','p','q'},'Before',1)
    
%% save fig
if  plotops.save_tif
    print(hfig, plotops.save_fig_filename, '-dpng', ['-r', num2str(plotops.output_resolution)]); % save the image as raster
%     exportgraphics(hfig,[plotops.save_fig_filename,'.tif'], 'Resolution',plotops.output_resolution); % save as raster
end
if  plotops.save_pdf
    exportgraphics(hfig,[plotops.save_fig_filename,'.pdf'], 'BackgroundColor','none', 'ContentType','vector'); % save as vector
end

