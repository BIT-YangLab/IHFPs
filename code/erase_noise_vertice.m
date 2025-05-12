function [ind_label, flg] = erase_noise_vertice(ind_label, mesh)
% DESCRIPTION:
%   Remove noise vertice in surface
%
% USAGE: 
%   [ind_label, flg] = erase_noise_vertice(ind_label, 'fs_LR_32k')
%
% Inputs:   ind_label: label of atlas 
%
%           mesh: fs_LR_32k or fsaverage
%
% Outputs:  ind_label: processed atlas
%
%           flg: 1 if there is parcel missing after removing noise, 0 otherwise
%
% Written by Jinlong Li, Guoyuan Yang
    %% prepare
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

    flg = 0;
    lh_vn = lh_targ_mesh.vertexNbors ;
    rh_vn = rh_targ_mesh.vertexNbors ;
    rh_vn(rh_vn ~= 0) = rh_vn(rh_vn ~= 0) + lh_end_idx;
    neighbors = [lh_vn rh_vn];

    vertex_n_ths = 10;
    new_label = find_local_consistent_region(ind_label, mesh);
    label_list = unique(new_label);
    label_list(label_list < 1) = [];

    label_cnt = max(label_list) ;
    for i = 1:length(label_list)
        index = find(new_label == label_list(i));
        net_tmp = ind_label(index(1));
        index_1 = find(ind_label == net_tmp);

        if length(index) < vertex_n_ths || ( isempty(setdiff(index_1, index)) && length(index) < length(setdiff(index_1, index)))
            ind_label(index) = label_cnt + i;
        end
    end

    label_list = unique(ind_label);
    label_list(label_list <= label_cnt) = [];
    
    for i = 1:length(label_list)
        par_index = find(ind_label == label_list(i));
        isvisit = zeros(size(par_index));
        while nnz(isvisit) < length(par_index)
            k_cnt = 0;
            for j = 1:length(par_index)
                if isvisit(j) == 1
                    continue;
                end
                nei_index = neighbors(:, par_index(j));
                nei_index(nei_index < 1) = [];
                nei_label = ind_label(nei_index);
                nei_label(nei_label == label_list(i)) = [];
                nei_label(nei_label < 1) = [];
                if length(nei_label) > 0
                    ind_label(par_index(j)) = mode(nei_label);
                    isvisit(j) = 1;
                else
                    k_cnt = k_cnt + 1;
                end
            end
            if k_cnt == length(par_index)
                flg = 1;
                return
            end
        end
    end

end