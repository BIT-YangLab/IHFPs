function roi_center_ind = compute_cortical_centroid_pos(labels, label_list, vertices_pos)
% DESCRIPTION:
%   Calculate the centroid of ROI on the cortex
%
% USAGE: 
%   label_list = unique(labels); label_list(label_list < 1) = [];
%   roi_center_ind = compute_cortical_centroid_pos(labels, label_list, vertices_pos);
%
% Inputs:   labels: labels of parcellation map  N vertex
%
%           label_list: list of labels
%        
%           vertices_pos: [x, y, z] position of vertex, N * 3 
%
% Outputs:  roi_center_ind:  index of centroid of ROI
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

roi_center_ind = zeros(length(label_list), 1);
for i = 1: length(label_list)

    roi_vertices_ind = find(labels==label_list(i));
    if isempty(roi_vertices_ind)
        continue
    end
    roi_vertices_pos = vertices_pos(labels==label_list(i), :);

    dist_x = repmat(roi_vertices_pos(:, 1), 1, size(roi_vertices_pos, 1)) - repmat(roi_vertices_pos(:, 1), 1, size(roi_vertices_pos, 1))';
    dist_y = repmat(roi_vertices_pos(:, 2), 1, size(roi_vertices_pos, 1)) - repmat(roi_vertices_pos(:, 2), 1, size(roi_vertices_pos, 1))';
    dist_z = repmat(roi_vertices_pos(:, 3), 1, size(roi_vertices_pos, 1)) - repmat(roi_vertices_pos(:, 3), 1, size(roi_vertices_pos, 1))';

    dist_W = sqrt(dist_x.^2 + dist_y.^2 + dist_z.^2);

    [~, min_sumdist_pos] = min(sum(dist_W, 1));
    
    roi_center_ind(i) = roi_vertices_ind(min_sumdist_pos);
end



end