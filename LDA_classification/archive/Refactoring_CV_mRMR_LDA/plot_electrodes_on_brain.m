function surf_plot_data_table = plot_electrodes_on_brain(edb_nonan_subset, class_label, ops)

% ops: optional struct containing fields:
%       colormap: either predefined matlab colormap or custom 256x3 RGB colormap...
%           ..... see https://www.mathworks.com/help/matlab/ref/colormap.html
%       clusters_as_colors: if true, assign colors based on cluster, using electrodes.global_clust_num... 'color_values' will be ignored
%          ... cluster colors will be consistent within a colormap and across alignment conditions,
%          ... based on Kuzdeba thesis table 1 start time
%       base_surfaces: list of surfaces to plot points on top of (cell array of strings); defaults to the 4 surf_show defaults
%       tempfile: .txt file for saving mni coordinates; this file will be deleted after being loaded by surf_show
%       img_savename: if not empty, the brain image will be saved as a .png with this filename
%       output_resolution: if img_savename is not empty, the saved image will be saved as this resolution in dots per inch (default 300)
%       force_close: if true, immediately close the surf_show window with 'close all force' (useful running this func in a batch to save figures)
%
%%% updated 2021/12/ by Liam Jackson

addpath /project/busplab/software/ecog/util
addpath /project/busplab/software/spm12
addpath /project/busplab/software/conn
addpath /project/busplab/software/display	
addpath /project/busplab/software/display/surf

lh_tempfile_default_name = ['lh_mnitemp_to_delete_', datestr(datetime), '.txt']; % use new filename each time
rh_tempfile_default_name = ['rh_mnitemp_to_delete_', datestr(datetime), '.txt']; % use new filename each time

vardefault('ops',struct);
field_default('ops','clusters_as_colors', 0)
field_default('ops','subjects_as_shapes', 0)
field_default('ops','suppress_brain_surfs', 0)
field_default('ops','base_surfaces', {'lh.pial.smoothed.surf', 'rh.pial.smoothed.surf','rh.subcortical.surf','lh.subcortical.surf'}); 
field_default('ops','lh_tempfile', ['/projectnb2/busplab/surfshowfiles/', lh_tempfile_default_name])
field_default('ops','rh_tempfile', ['/projectnb2/busplab/surfshowfiles/', rh_tempfile_default_name])
field_default('ops','img_savename',[])
field_default('ops','output_resolution',300); 
field_default('ops','selected_hemi','both');
field_default('ops','viewpoint','left'); % anterior posterior left right superior inferior
field_default('ops','do_mosaic', 0);
field_default('ops','force_close',0); 
field_default('ops','top_k_pct',[]); 

clusterkey_file =  '/project/busplab/software/ecog/data/clusterkey'; 
load(clusterkey_file, 'clusterkey');

surfshow_resolution = 1; % 1= lowres, 2 = highres

num_electrodes = height(edb_nonan_subset);

% if toggled, map color to cluster number
if ops.clusters_as_colors
    close all % trying to make the legend with other figs open can cause problems
    color_values = edb_nonan_subset.global_clust_num; % get cluster numbers of selected electrodes
    mincolval = 1; % fix min and max vals to encompass all 8 cluster numbers
    maxcolval = 8; % fix min and max vals to encompass all 8 cluster numbers
    % make colormap for clusters_as_colors that shows each of the 8 clusters distinctly
    clusts_cmap = [ 1.0     0.5     0.0 ;...
                    0.9     0.9     0.0 ;...
                    0.1     1.0     0.0 ;...
                    0.0     0.7     1.0 ;...
                    0.0     0.0     1.0 ;...
                    0.6     0.5     1.0 ;...
                    0.8     0.0     0.8 ;...
                    1.0     0.0     0.0 ];
    field_default('ops','colormap', clusts_cmap);
    color_values = clusts_cmap(color_values, :);
else
    red_rgb = [1 0 0]; 
    color_values = repmat(red_rgb, num_electrodes, 1);
    clusts_cmap = jet(8);
    field_default('ops','colormap', clusts_cmap);
end

collapse_this = 0;
while collapse_this
% make sure there is exactly 1 color value for every electrode
% assert(num_electrodes == size(color_values, 1), 'electrode_table_indices must be same length as color_values')

% map color vals to RGB values in chosen color map
% vardefault('mincolval', min(color_values)); % set value range if not already fixed by clusters_as_colors
% vardefault('maxcolval', max(color_values)); % set value range if not already fixed by clusters_as_colors
% normvals = (color_values - mincolval) ./ (maxcolval-mincolval); % fit vals to range 0-1... baselined vals minus range
% rgbinds = 1 + round(normvals * [size(clusts_cmap,1)-1]); % match normed vals to indices in colormap
% rgbvals = ops.colormap(rgbinds,:); 
end

shapes_options = {  's',...    % - sphere (size --> side length, in mm) 'u',...    % - cube in MNI space (size --> diameter, in mm) || 'b',...    % - cube in image space (size --> diameter, in voxels %                     'c',...    % - crosshair in MNI space (size --> length of each axis, in mm || 'r',...    % - crosshair in image space (size --> length of each axis, in voxels
                    'v',...    % - plumb bob
                    'z',...    % - tetrahedron (triangular pyramid) % 't',...    % - tetrahedron (rotated pyramid)
                    'a',...    % - asterisk
                    'i'};      % - icosahedron
shapes_options_verbose = {  'sphere',... % 'cube',... % 'cross',... 
                            'diamond',... % actually plumb bob
                            'tetrahedron',... 
                            'asterisk',...
                            'hexagon'}; % actually icosahedron
                
sub_nums_vector = edb_nonan_subset.subject;
shapes = char(num_electrodes, 1);
shapes_verbose = cell(num_electrodes, 1);

if ops.subjects_as_shapes                
    shapes = char(edb_nonan_subset.surf_shape);
    for s_idx = 1:length(shapes)
        shapes_verbose(s_idx) = shapes_options_verbose(strcmp(shapes(s_idx), shapes_options));
    end
else
    shapes(:) = shapes_options{1};
end

lefthemi_vector = isnan(edb_nonan_subset.right_hemi);
righthemi_vector = ~isnan(edb_nonan_subset.right_hemi);
alignment_vector = edb_nonan_subset.type;
global_clust_num_vector = edb_nonan_subset.global_clust_num; 
surf_plot_data_table = table(alignment_vector, sub_nums_vector, shapes, shapes_verbose, global_clust_num_vector,...
    'VariableNames', {'type', 'sub_num', 'shape', 'shape_long', 'global_clust_num'});
% clust_name_vector = cell(height(surf_plot_data_table),1);
% clust_name_vector(:) = clusterkey.name(clusterkey.global_clust_num == global_clust_num_vector(1));
% surf_plot_data_table.clust_name = clust_name_vector; 

% save mni coordinates to be plotted as .txt file
mni_coords = edb_nonan_subset.mni_coord; 

% writematrix([mni_coords, zeros(size(electrode_table_indices)), shapes], ops.tempfile);
mni_coords = mat2cell(mni_coords, ones(size(mni_coords, 1),1), ones(1,size(mni_coords, 2)));
mni_coords(:,4) = mat2cell(ones(size(num_electrodes)), ones(size(num_electrodes, 1), 1), ones(1, size(num_electrodes, 2)));
mni_coords(:,5) = mat2cell(shapes, ones(size(shapes, 1), 1), ones(1, size(shapes, 2)));
writecell(mni_coords(lefthemi_vector, :), ops.lh_tempfile);
writecell(mni_coords(righthemi_vector, :), ops.rh_tempfile);

% 
% if ~ops.suppress_brain_surfs
% plot electrode points with specified colors on brain surface
alpha = 1;

if strcmp(ops.viewpoint, 'medial')    
    switch ops.selected_hemi
        case 'both'
            error('Invalid selection. Medial viewpoint selected as well as both hemispheres');
        case 'left'
            view = 'right view';
        case 'right'
            view = 'left view';
    end
else
    view = [ops.viewpoint, ' view'];
end

if ops.do_mosaic
    horz_off = 150;
    left_shift = [0, horz_off, 0]; 
    right_shift = [0, -horz_off, 0];
else
    left_shift = [0 0 0]; 
    right_shift = [0 0 0];
end

if strcmp(ops.selected_hemi, 'both')
surf_show(  'SURFACE_ADD', ops.base_surfaces{1}, 'SURFACE_PROPERTIES', 'offset', right_shift, 'transparency', alpha,... % LH pial base surface
            'SURFACE_ADD', ops.base_surfaces{4}, 'SURFACE_PROPERTIES', 'offset', right_shift, 'transparency', alpha,... % LH subcortical base surface
            'SURFACE_ADD', ops.base_surfaces{2}, 'SURFACE_PROPERTIES', 'offset', left_shift, 'transparency', alpha,... % RH pial base surface
            'SURFACE_ADD', ops.base_surfaces{3}, 'SURFACE_PROPERTIES', 'offset', left_shift, 'transparency', alpha,... % RH subcortical base surface
            'SURFACE_ADD', ops.rh_tempfile, 'SURFACE_PROPERTIES','offset', left_shift, 'color', color_values(righthemi_vector,:),...
            'SURFACE_ADD', ops.lh_tempfile, 'SURFACE_PROPERTIES','offset', right_shift, 'color', color_values(lefthemi_vector,:), 'VIEW', view, 'RESOLUTION', surfshow_resolution); %, 'PRINT', 'view', 'mosaic'); % colored electrode locations

elseif strcmp(ops.selected_hemi, 'left') && ~ops.do_mosaic
surf_show(  'SURFACE_ADD', ops.base_surfaces{1}, 'SURFACE_PROPERTIES', 'offset', right_shift, 'transparency', alpha,... % LH pial base surface
            'SURFACE_ADD', ops.base_surfaces{4}, 'SURFACE_PROPERTIES', 'offset', right_shift, 'transparency', alpha,... % LH subcortical base surface
            'SURFACE_ADD', ops.lh_tempfile, 'SURFACE_PROPERTIES','offset', right_shift, 'color', color_values(lefthemi_vector,:), 'VIEW', view, 'RESOLUTION', surfshow_resolution); %, 'PRINT', 'view', 'mosaic'); % colored electrode locations

elseif strcmp(ops.selected_hemi, 'right') && ~ops.do_mosaic
surf_show(  'SURFACE_ADD', ops.base_surfaces{2}, 'SURFACE_PROPERTIES', 'offset', left_shift, 'transparency', alpha,... % RH pial base surface
            'SURFACE_ADD', ops.base_surfaces{3}, 'SURFACE_PROPERTIES', 'offset', left_shift, 'transparency', alpha,... % RH subcortical base surface
            'SURFACE_ADD', ops.rh_tempfile, 'SURFACE_PROPERTIES','offset', left_shift, 'color', color_values(righthemi_vector,:), 'VIEW', view, 'RESOLUTION', surfshow_resolution); %, 'PRINT', 'view', 'mosaic'); % colored electrode locations

else
    error('Something is wrong');
end
delete(findobj('Type','Line')) % delete dotted outlines around patches
% else
%     surf_show(  'SURFACE_ADD', ops.rh_tempfile, 'SURFACE_PROPERTIES','color', color_values(righthemi_vector,:), ...
%                 'SURFACE_ADD', ops.lh_tempfile, 'SURFACE_PROPERTIES','color', color_values(lefthemi_vector,:), 'RESOLUTION', surfshow_resolution); %, 'PRINT', 'view', 'mosaic'); % colored electrode locations
%     delete(findobj('Type','Line')) % delete dotted outlines around patches
% end

add_legend(ops, clusterkey, surf_plot_data_table)
add_title(ops, surf_plot_data_table, class_label)

% if specified, save 
if ~isempty(ops.img_savename)
    surfshow_ax= gca;
    set(surfshow_ax,'Position',[.15 .15 .7 .7])
    tempfig= figure('Visible','off'); % Invisible figure
    temp_ax = copyobj(surfshow_ax,tempfig); % Copy the appropriate axes
    add_legend(ops, clusterkey, surf_plot_data_table)
    add_title(ops, surf_plot_data_table, class_label)
    print(tempfig, ops.img_savename, '-dpng', ['-r', num2str(ops.output_resolution)]); % save the image file
    close(tempfig)
end
    
delete(ops.lh_tempfile) % delete temporary mni coordinates file
delete(ops.rh_tempfile) % delete temporary mni coordinates file
if ops.force_close
    close all force % close the surf_show window without requiring manual approval
end

end


function add_legend(ops, clusterkey, t)
if ops.clusters_as_colors
    add_clusters_legend(ops, clusterkey)
end
if ops.subjects_as_shapes
    add_shapes_legend(t)
end

function add_clusters_legend(ops, clusterkey)
    allpatches = table(findobj('Type','Patch'), 'VariableNames',{'patch'});
    npatches = height(allpatches);
    allpatches.global_clust_num = NaN(npatches,1); 
    % in surf show, skip first 2 and last 2 Patch objects, which are something other than the points we added
%     for ipatch = 3:npatches-2
    for ipatch = 1:npatches
        thiscol = allpatches.patch(ipatch).FaceVertexCData(1,:); % color of this patch 
        thisclust = find(all(thiscol == ops.colormap,2)); 
        if ~isempty(thisclust)
            allpatches.global_clust_num(ipatch) = thisclust; % find which clustnum produced this patch color
        end
    end
    [clustnums, inds4leg] = unique(allpatches.global_clust_num);
    inds4leg = inds4leg(~isnan(clustnums)); % exclude patches which are not points we plotted
    clustnums = clustnums(~isnan(clustnums)); 
    nclusts = length(clustnums);
    leglabels = cell(size(clustnums)); 
    for iclust = 1:nclusts % look up the cluster names from clusterkey
        clusterkey_matchrow = find(clusterkey.global_clust_num==clustnums(iclust),1);
        leglabels{iclust} = clusterkey.name{clusterkey_matchrow}; 
    end
    legend(allpatches.patch(inds4leg),leglabels, 'Position', [.85 .16 .05 .05])
end

function add_shapes_legend(t)
    t = unique(t(:,{'sub_num', 'shape_long'}));

    shapes_long = t.shape_long;
    shapes_long_proper = cell(length(shapes_long), 1);
    for i = 1:length(shapes_long)
        shapes_long_proper(i) = {plaintext(shapes_long{i})};
    end

    sub_nums_str = string(t.sub_num);
    ss = cell(length(sub_nums_str), 1);
    ss(:) = {'S'};
    colons = {': '};
    manleg = strcat(ss, sub_nums_str, colons, shapes_long_proper);

    fig = gcf;
    fig.Units = 'normalized';
    ax = gca; 

    left = ax.Legend.Position(1); 
    bottom = ax.Legend.Position(2); 
    w = ax.Legend.Position(3); 
    h = ax.Legend.Position(4); 

    annotation('textbox', [0.075, 0.16, .05, .05], 'String', manleg, 'FitBoxToText', 'on')

end
end

function add_title(ops, t, class_label)
type = t.type(1); 
switch type
    case 2
        alignment = 'Stimulus';
    case 1
        alignment = 'Onset'; 
end
% clust_name = t.clust_name{1}; 
if ~isempty(ops.top_k_pct)
    title_str = sprintf('Top %2.0f%% %s Data, %s', ops.top_k_pct * 100, alignment, plaintext(class_label));
else
    title_str = sprintf('%s Data, %s', alignment, plaintext(class_label));
end
title({title_str, sprintf('%s Hemi, %s view', plaintext(ops.selected_hemi), plaintext(ops.viewpoint))}, 'FontSize', 18)

end