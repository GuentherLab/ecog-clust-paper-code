function job = execute_multinode_feat_importance(params)

number_of_nodes = params.number_of_nodes;

% Note: temp_params_filename is needed because if you just pass conn the
% params object, the 10gb of local user space fills and the pipeline
% shuts down. To avoid this, saving a temp filename and pass that to the
% conn function in a struct. cv_mRMR_LDA_SCC_node_pipeline then loads the
% params obj from the temporary filename 
temp_params_filename = fullfile(params.grouping_path, params.current_group_value, 'params.mat'); 
save(temp_params_filename, 'params', '-v7.3');

addpath /project/busplab/software/conn 
addpath /project/busplab/software/spm12

temp_struct = struct(); 
temp_struct.temp_params_filename = temp_params_filename; 

for node_number = 1:number_of_nodes 
    temp_struct.node_number = node_number; 
%     multinode_feat_importance(temp_struct)
    job(node_number) = conn('submit', 'multinode_feat_importance', temp_struct);
end

end

