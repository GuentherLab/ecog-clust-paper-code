 %%% plot top encoding elelctrodes for subject (362) with electrodes in both hemispheres
% for supp fig

clear

setelc_ops.sub = 362; 
setelc_ops.proportion_top_elcs = 0.6; 
setelc_ops.surf_shape = 's';
setelc_ops.colormap = copper; 

surf_ops.scale = 0.8; % electrode marker size
surf_ops.force_close = 0; % if true, close the surf_show window after plotting/saving image
surf_ops.crop_borders_top_bot_left_right = [230,1100,240,1500]; % if not empty, crop saved image borders at these pixel values

%%%%%%%% options for normalizing encoding strengths
setelc_ops.strength_norming = 'within_phonunit';
% setelc_ops.strength_norming = 'within_phonunit_baseline_zero';
% setelc_ops.strength_norming = 'across_phonunit'; % max value in colorbar will also be set to maximum strength across all phonunits

set_paths()
if ~exist('elc','var')
    load(TOP_ELC_DATA_FILE)
    elc.Properties.VariableNames = strrep(elc.Properties.VariableNames,'word','syl');
end
setelc_ops.output_dir = [DATA_DIR filesep 'figs' filesep 'top_elcs_single_sub'];

surf_ops.lh_tempfile = [DATA_DIR filesep 'lh_mnitemp_to_delete_', strrep(datestr(datetime),':','-'), '.txt']; % use new filename each time
surf_ops.rh_tempfile = [DATA_DIR filesep 'rh_mnitemp_to_delete_', strrep(datestr(datetime),':','-'), '.txt']; % use new filename each time
surf_ops.surfshowDir = 'C:/docs/code/matlab/display/surf';
surf_ops.clusterkey_file = [SOURCEDATA_DIR filesep 'clusterkey.mat'];
% surf_ops.viewpoint = 'lateral';
surf_ops.add_legend = 0; 
% surf_ops.colormap = [0 0 0];

mkdir(setelc_ops.output_dir );



%%
setelc_ops.phonunit = 'cons';
surf_ops.selected_hemi = 'left'; 
surf_ops.viewpoint = 'left';
[elc_temp_cons, surf_ops] = set_elcs(elc,setelc_ops,surf_ops); 
surftab = plot_electrodes_on_brain_AM2(elc_temp_cons,surf_ops);


%%
setelc_ops.phonunit = 'cons';
surf_ops.selected_hemi = 'right'; 
surf_ops.viewpoint = 'left';
[elc_temp_cons, surf_ops] = set_elcs(elc,setelc_ops,surf_ops); 
surftab = plot_electrodes_on_brain_AM2(elc_temp_cons,surf_ops);

%%
setelc_ops.phonunit = 'cons';
surf_ops.selected_hemi = 'right'; 
surf_ops.viewpoint = 'right';
[elc_temp_cons, surf_ops] = set_elcs(elc,setelc_ops,surf_ops); 
surftab = plot_electrodes_on_brain_AM2(elc_temp_cons,surf_ops);

%%
setelc_ops.phonunit = 'vowel';
surf_ops.selected_hemi = 'left'; 
surf_ops.viewpoint = 'left';
[elc_temp_vowel, surf_ops] = set_elcs(elc,setelc_ops,surf_ops); 
surftab = plot_electrodes_on_brain_AM2(elc_temp_vowel,surf_ops);

%%
setelc_ops.phonunit = 'vowel';
surf_ops.selected_hemi = 'right'; 
surf_ops.viewpoint = 'left';
[elc_temp_vowel, surf_ops] = set_elcs(elc,setelc_ops,surf_ops); 
surftab = plot_electrodes_on_brain_AM2(elc_temp_vowel,surf_ops);

%%
setelc_ops.phonunit = 'vowel';
surf_ops.selected_hemi = 'right'; 
surf_ops.viewpoint = 'right';
[elc_temp_vowel, surf_ops] = set_elcs(elc,setelc_ops,surf_ops); 
surftab = plot_electrodes_on_brain_AM2(elc_temp_vowel,surf_ops);

%%
setelc_ops.phonunit = 'syl';
surf_ops.selected_hemi = 'left'; 
surf_ops.viewpoint = 'left';
[elc_temp_syl, surf_ops] = set_elcs(elc,setelc_ops,surf_ops); 
surftab = plot_electrodes_on_brain_AM2(elc_temp_syl,surf_ops);

%%
setelc_ops.phonunit = 'syl';
surf_ops.selected_hemi = 'right'; 
surf_ops.viewpoint = 'left';
[elc_temp_syl, surf_ops] = set_elcs(elc,setelc_ops,surf_ops); 
surftab = plot_electrodes_on_brain_AM2(elc_temp_syl,surf_ops);

%%
setelc_ops.phonunit = 'syl';
surf_ops.selected_hemi = 'right'; 
surf_ops.viewpoint = 'right';
[elc_temp_syl, surf_ops] = set_elcs(elc,setelc_ops,surf_ops); 
surftab = plot_electrodes_on_brain_AM2(elc_temp_syl,surf_ops);



%% create colorbars
% save each colorbar fig as .svg then edit in inkscape
nticks = 5;

close all force

for this_phonunit = {'cons','vowel','syl'}
    % close all 
    hfig = figure;
    hax = gca; 
    hax.Visible = 'off';
    hcb = colorbar;
    colormap(setelc_ops.colormap)
    eval(['thistab = elc_temp_', this_phonunit{:}, ';']) % copy topelc table for this phonunit

    phonunit_strengths = thistab{:, [this_phonunit{:} '_accuracy_change_wo_electrode']};
    max_phonunit_strength = max(phonunit_strengths)
    
    %%% label ticks in percent accuracy improvement
    hcb.Ticks = linspace(0,1,nticks);
    hcb.TickLabels = cellstr(num2str([round(1000*linspace(0,max_phonunit_strength,nticks))/10]'));
    hcb.FontName = 'Arial';


end





%% statistical test on R vs. L hemisphere effects
setelc_ops_for_stats = setelc_ops;

setelc_ops_for_stats.proportion_top_elcs = 0.6; 
surf_ops.selected_hemi = 'both'; 

% setelc_ops_for_stats.outcome_var = 'acc_change'; 
setelc_ops_for_stats.outcome_var = 'acc_change_normed'; 
% setelc_ops_for_stats.outcome_var = 'rank_within_phonunit'; 

setelc_ops_for_stats.strength_norming = 'within_phonunit';
% setelc_ops_for_stats.strength_norming = 'within_phonunit_baseline_zero';
% setelc_ops_for_stats.strength_norming = 'across_phonunit'; % max value in colorbar will also be set to maximum strength across all phonunits

%
setelc_ops_for_stats.phonunit = 'cons'; [elc_temp_cons, surf_ops] = set_elcs(elc,setelc_ops_for_stats,surf_ops); 
setelc_ops_for_stats.phonunit = 'vowel'; [elc_temp_vowel, surf_ops] = set_elcs(elc,setelc_ops_for_stats,surf_ops); 
setelc_ops_for_stats.phonunit = 'syl'; [elc_temp_syl, surf_ops] = set_elcs(elc,setelc_ops_for_stats,surf_ops); 

topelc_1sub = table;
for this_phonunit = {'cons', 'vowel', 'syl'}
    eval(['thistab = elc_temp_', this_phonunit{:}, ';']) % copy topelc table for this phonunit
    topelc_1sub = [topelc_1sub; thistab(:,{'phonunit','acc_change','acc_change_normed','rank_within_phonunit','right_hemi','electrode'})];
end

topelc_1sub.hemi = cell(height(topelc_1sub),1);
topelc_1sub.hemi(isnan(topelc_1sub.right_hemi)) = {'l'};
topelc_1sub.hemi(~isnan(topelc_1sub.right_hemi)) = {'r'};
topelc_1sub.right_hemi = [];

% use accuracy change as the outcome variable; use phonunit and hemisphere as predictor
[p,tbl,stats,terms] = anovan(topelc_1sub{:,setelc_ops_for_stats.outcome_var }, {topelc_1sub.phonunit, topelc_1sub.hemi},'model',2,'varnames',{'phonunit','hemi'});

coeftable = table(stats.coeffnames,stats.coeffs,'VariableNames',{'coefname','coef'});




%%
function [elc_temp, surf_ops] = set_elcs(electrodes,setelc_ops,surf_ops)
    field_default('surf_ops','selected_hemi','both');

    if isfield(setelc_ops,'sub') && ~isempty(setelc_ops.sub)  % use specific sub
        elc_temp = electrodes(electrodes.subject==setelc_ops.sub,:);
    else     % use all subs
        elc_temp = electrodes;
    end

    max_encoding_strength_across_phonunits = max([elc_temp.cons_accuracy_change_wo_electrode;...
                                                 elc_temp.vowel_accuracy_change_wo_electrode;...
                                                 elc_temp.syl_accuracy_change_wo_electrode]); 

    rankvar = [setelc_ops.phonunit '_accuracy_change_wo_electrode'];
    [~, elc_encoding_rank] = sort(elc_temp{:,rankvar},'descend'); % 'descend' b/c highest value = best encoding
    elc_temp = elc_temp(elc_encoding_rank, :); % sort rows by encoding strength on chosen phon unit
    n_elcs_to_keep = round(height(elc_temp) * setelc_ops.proportion_top_elcs);
    elc_temp = elc_temp(1:n_elcs_to_keep,:); % discard all but top encoders

    % assign colors based on decoding strength
    switch setelc_ops.strength_norming
        case 'within_phonunit_baseline_zero'     % set zero to rankvar minimum rather than absolute zero
            elc_temp{:,[rankvar '_normed']} = [elc_temp{:,rankvar} - min(elc_temp{:,rankvar}) ] / range(elc_temp{:,rankvar}); 
        case 'within_phonunit'
             elc_temp{:,[rankvar '_normed']} = elc_temp{:,rankvar} / max(elc_temp{:,rankvar});
        case 'across_phonunit'% norm across phonunits
            elc_temp{:,[rankvar '_normed']} = elc_temp{:,rankvar} / max_encoding_strength_across_phonunits; 
    end
    elc_temp.colormap_ind = round(255 * elc_temp{:,[rankvar '_normed']} + 1); % index in 256-value colormap
    elc_temp.color = setelc_ops.colormap(elc_temp.colormap_ind, :);
    elc_temp = movevars(elc_temp,{'electrode',rankvar,'right_hemi',[rankvar '_normed'],'color'},'Before',1);

    % add vars for later stats tests
    elc_temp.phonunit = cell(height(elc_temp),1);
    elc_temp.phonunit(:) = {setelc_ops.phonunit};
    elc_temp.acc_change = elc_temp{:,rankvar}; 
    elc_temp.acc_change_normed = elc_temp{:,[rankvar '_normed']}; 
    elc_temp.rank_within_phonunit = [1:height(elc_temp)]';

    % if a hemisphere is selected, exclude electrodes from the other hemisphere....
    % ....... so they don't accidentally get plotted
    if   surf_ops.selected_hemi ~= "both"
        if surf_ops.selected_hemi == "left"
            elc_temp = elc_temp(isnan(elc_temp.right_hemi),:);
        elseif surf_ops.selected_hemi == "right"
            elc_temp = elc_temp(~isnan(elc_temp.right_hemi),:);
        end
   end

    surf_ops.color_values_override = elc_temp.color;

    field_default('surf_ops', 'viewpoint', surf_ops.selected_hemi);

    surf_ops.img_savename = [setelc_ops.output_dir filesep 'sub_' num2str(setelc_ops.sub) '_brain_top_'...
        num2str(round(100 * setelc_ops.proportion_top_elcs)) 'pct_' setelc_ops.phonunit '_'  surf_ops.selected_hemi,...
        '_view_', surf_ops.viewpoint];
end