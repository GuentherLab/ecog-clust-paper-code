  %%% create plots representing the two patterns of hypothesized syllable and phoneme representation
 %
% updated by AM 2022/6/15
clear 

%% set timepoints
 ops.t_onset = 0.7; % voice onset time in sec
  ops.speech_dur = 1; % syllable production duration
 ops.post_speech_dur = .5; % time between speech offset and trial end
 
 
 %%  specify 6 values representing encoding strength:
 % 1. trial start
 % 2. speech onset 1
 % 3. speech onset 2
 % 4. speech offset 1
 % 5. speech offset 2
 % 6. end of trial

  % ESE = Early Syllable Encoding hypothesis
 ese_syl_tcourse = [1, 1, 0.3, 0.3, 0, 0];
ese_phon_tcourse = [0.1, 0.1, 0.8, 0.8, 0, 0]; 

 % EPE = Early Phoneme Encoding hypothesis
 epe_syl_tcourse = [0.1, 0.1, 0.3, 0.3, 0, 0];
epe_phon_tcourse = [1, 1, 1, 0, 0, 0];  
 
%% plot params
 ops.ylimits = [0 1.2]; 
% ops.xlimits = [0, 7];
% ops.leg_location = 'east';
ops.xlines_width = 1;
    ops.xline_style = '--';
    ops.xline_color =  [0.15 0.15 0.15];
% ops.annot_pos_xy =  [0.70, 0.78]; % position of pval annotation
ops.line_color = [0.2 0.2 0.2]; 
ops.line_width = 3;
ops.axes_line_width =  1;
ops.axis_font_size =  13;
ops.axes_numbers_bold =  'bold';
ops.font =  'Arial';
ops.fig_width_length =  [900 600];
ops.box_on_off = 'off';
ops.yticks = [];
ops.xticks = [];
% ops.xticklength = [0 0];
% ops.x_tick_angle =  0; % in degrees; 0=straight; positive rotates counter-clockwise
ops.background_color =  [1 1 1]; 

 
 
 %% make plots
 t_offset = ops.t_onset +  ops.speech_dur;
 t_trial_end = t_offset + ops.post_speech_dur;
 ops.xvals = [0, ops.t_onset, ops.t_onset, t_offset, t_offset, t_trial_end];
 

 close all
 
%  hfig = figure('Windowstyle', 'undocked');
     hfig = figure('Windowstyle', 'docked');
 hfig.Color = ops.background_color;
 set(hfig,'Renderer', 'painters', 'Position', [50, 50, ops.fig_width_length(1), ops.fig_width_length(2)])

 
subplot(2,2,1)
hplot = plot(ops.xvals, ese_phon_tcourse);
hax = gca; 
hax = plotformat(hax, ops);
ylabel({'Phoneme', 'encoding strength'})

subplot(2,2,2)
hplot = plot(ops.xvals, ese_syl_tcourse);
hax = gca; 
hax = plotformat(hax, ops);
ylabel({'Syllable', 'encoding strength'})

subplot(2,2,3)
hplot = plot(ops.xvals, epe_phon_tcourse);
hax = gca; 
hax = plotformat(hax, ops);
ylabel({'Phoneme', 'encoding strength'})

subplot(2,2,4)
hplot = plot(ops.xvals, epe_syl_tcourse);
hax = gca; 
hax = plotformat(hax, ops);
ylabel({'Syllable', 'encoding strength'})



%% plot ops function
function haxout = plotformat(hax, ops)
      ylim(ops.ylimits)
    box(ops.box_on_off)


    hax.XTick = ops.xticks;
%     set(hax,'XTickLabels', clustlist)
%     set(hax,'XColor',[0 0 0]);
%     set(hax,'YColor',[0 0 0]);
    hax.ColorOrder = ops.line_color;
    hax.Children(1).LineWidth = ops.line_width;
    set(hax,'linewidth', ops.axes_line_width)
    set(hax,'FontSize', ops.axis_font_size)
    set(hax,'FontWeight', ops.axes_numbers_bold)
    set(hax,'FontName', ops.font)
    set(hax,'YTick',ops.yticks)
%     xlim([ops.xlimits]);
%     h=gca; h.XAxis.TickLength = ops.xticklength;
%     hylabel = ylabel({'Proportion of top', 'Word or Consonant electrodes'});
%         hylabel.Position = [-0.7 0.5 0.5];
    xlabel('Time')
    hxline = xline(ops.xvals(2));
        hxline.LineWidth = ops.xlines_width;
        hxline.LineStyle = ops.xline_style;
        hxline.Color = ops.xline_color;
    hxline = xline(ops.xvals(4));
        hxline.LineWidth = ops.xlines_width;
        hxline.LineStyle = ops.xline_style;
        hxline.Color = ops.xline_color;
    haxout = hax;
  end
    
    
  