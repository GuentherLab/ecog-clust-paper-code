%     fprintf(1, '    s - sphere (size --> side length, in mm)\n');
%     fprintf(1, '    c - crosshair in MNI space (size --> length of each axis, in mm\n');
%     fprintf(1, '    r - crosshair in image space (size --> length of each axis, in voxels\n');
%     fprintf(1, '    u - cube in MNI space (size --> diameter, in mm)\n');
%     fprintf(1, '    b - cube in image space (size --> diameter, in voxels)\n');

ops = struct(); 
ops.force_close = 0; 
ops.clusters_as_colors = 1; 
ops.subjects_as_shapes = 1; 

electrode_table_indices = out.global_elec_idx;
color_values = out.colors; 
alpha = out.alpha; 
shapes = out.shapes; 
shapes = shapes(1:10);
shapes(8) = 'c';
shapes(9) = 't';
shapes(10) = 'i';
shapes(11) = 't';
shapes(12) = 'i';
plot_greenlee_electrodes_on_brain(electrode_table_indices, color_values, alpha, shapes, ops)
%%
close all; 
f = 1; 
v = 2; 
j = length(a{f}(:,1)); 
k = length(a{v}(:,1)); 
% scatter3(a{f}(1:j, 1), a{f}(1:j, 2), a{f}(1:j, 3))
hold on ; 
surf(a{v}(1:k, 1), a{v}(1:k, 2), a{v}(1:k, 3))

%%
