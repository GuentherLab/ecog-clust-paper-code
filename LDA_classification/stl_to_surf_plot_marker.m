% This adds data from an STL file to the VTK_ReferenceShapes file. 
% Adds shape options to surf plots, must generate binary STL file in CAD 
% then add the option in /project/busplab/software/display/CreateVolume.m, 
% line 286-ish. Also has fxn to remove shapes from VTKfile.

stl_filename = 'triangular_pyramid.stl' ;
shape_name = 'tet2'; 

add_surf_shape(stl_filename, shape_name)

% Use this to verify the data was added correctly
vtk_mf = matfile('/project/busplab/software/display/VTK_ReferenceShapes.mat');

patches = vtk_mf.patches
shapes = vtk_mf.shapes

% remove_surf_shape('ast');

function add_surf_shape(stl_filename, shape_name)

stl = stlread(stl_filename); 

stl_patch = trimesh(stl);  
stl_faces = stl_patch.Faces;
stl_vertices = stl_patch.Vertices;

% you can use this section to scale any dimension individually

% x_scale = 2; 
% y_scale = 2; 
% z_scale = 2; 
global_scale = .75; 
% 
% stl_vertices(:,1) = stl_vertices(:,1) .* x_scale; % x-scaling
% stl_vertices(:,2) = stl_vertices(:,2) .* y_scale; % y-scaling
% stl_vertices(:,3) = stl_vertices(:,3) .* z_scale; % z-scaling
stl_vertices = stl_vertices .* global_scale; % global-scaling

close all;

vtk_mf = matfile('/project/busplab/software/display/VTK_ReferenceShapes.mat', 'Writable', true);

patches = vtk_mf.patches;
shapes = vtk_mf.shapes;

if ~ismember(shape_name, shapes)
    current_num_shapes = length(shapes); 
    shapes{current_num_shapes + 1} = shape_name; 

    current_num_patches = length(patches); 
    patches{current_num_patches + 1}.vertices = stl_vertices; 
    patches{current_num_patches + 1}.faces = stl_faces; 
else
    shape_idx = strcmp(shape_name, shapes); 
    patches{shape_idx}.vertices = stl_vertices; 
    patches{shape_idx}.faces = stl_faces;     
end
    
vtk_mf.patches = patches; 
vtk_mf.shapes = shapes;

end

function remove_surf_shape(shape_name)

vtk_mf = matfile('/project/busplab/software/display/VTK_ReferenceShapes.mat', 'Writable', true);

patches = vtk_mf.patches;
shapes = vtk_mf.shapes;

if ~ismember(shape_name, shapes)
    fprintf('Shape %s not found, nothing modified.\n', shape_name); 
else
    shape_idx = strcmp(shape_name, shapes); 
    patches(shape_idx) = [];
    shapes(shape_idx) = [];
    fprintf('Shape %s removed successfully.\n', shape_name); 
end
    
vtk_mf.patches = patches; 
vtk_mf.shapes = shapes;

end
