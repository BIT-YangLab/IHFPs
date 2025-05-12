function dist_W = compute_cortical_vertices_dist(vertices_i, vertices_j, vertices_pos)
% DESCRIPTION:
%   Calculate the Euclidean distance between each vertex
%
% USAGE: 
%   label_list = unique(labels); label_list(label_list < 1) = [];
%   dist_W = compute_cortical_vertices_dist(vertices_i, vertices_j, vertices_pos);
%
% Inputs:   vertices_i, vertices_j: verteces with index i,j 
%
%           vertex_pos: matrix of position  
%
% Outputs:  dist_W:  Euclidean distance between each vertex
%
% Written by Jinlong Li, Guoyuan Yang
    %% prepare
    if(~exist('vertices_pos', 'var'))
        mesh = 'fs_LR_32k';
        if contains(mesh, 'fs_LR')
            lh_targ_mesh = CBIG_read_fslr_surface('lh', mesh,'inflated','medialwall.annot');
            rh_targ_mesh = CBIG_read_fslr_surface('rh', mesh,'inflated','medialwall.annot');
            no_medial_wall_index = [find(lh_targ_mesh.MARS_label == 2); find(rh_targ_mesh.MARS_label == 2) + length(lh_targ_mesh.MARS_label)];
            lh_end_idx = 32492;
        elseif contains(mesh, 'fsaverage') 
            lh_targ_mesh = CBIG_ReadNCAvgMesh('lh', mesh, 'inflated', 'cortex');
            rh_targ_mesh = CBIG_ReadNCAvgMesh('rh', mesh, 'inflated', 'cortex');
            
            no_medial_wall_index = [find(lh_targ_mesh.MARS_label == 2) find(rh_targ_mesh.MARS_label == 2) + length(lh_targ_mesh.MARS_label)]';
            lh_end_idx = 40962;
        end
    
    
        vertices_pos = [lh_targ_mesh.vertices'; rh_targ_mesh.vertices'];
    end

    roi_vertices_pos_i = vertices_pos(vertices_i, :);
    roi_vertices_pos_j = vertices_pos(vertices_j, :);

    dist_x = repmat(roi_vertices_pos_i(:, 1), 1, size(roi_vertices_pos_j, 1)) - repmat(roi_vertices_pos_j(:, 1), 1, size(roi_vertices_pos_i, 1))';
    dist_y = repmat(roi_vertices_pos_i(:, 2), 1, size(roi_vertices_pos_j, 1)) - repmat(roi_vertices_pos_j(:, 2), 1, size(roi_vertices_pos_i, 1))';
    dist_z = repmat(roi_vertices_pos_i(:, 3), 1, size(roi_vertices_pos_j, 1)) - repmat(roi_vertices_pos_j(:, 3), 1, size(roi_vertices_pos_i, 1))';

    dist_W = sqrt(dist_x.^2 + dist_y.^2 + dist_z.^2);

end