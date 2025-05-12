function new_label_t = network_striped(group_label, network_match)
% DESCRIPTION:
%   Visualization of stripped network
%
% USAGE: 
%   new_label_t = network_striped(group_label, network_match)
%
% Inputs:   group_label: group atlas
%
%           network_match: 17x1 cell, parcel index corresponding to each network
%
% Outputs:  new_label_t: processed labels for visualization
% Written by Jinlong Li, Guoyuan Yang

lh_targ_mesh = CBIG_read_fslr_surface('lh', 'fs_LR_32k','inflated','medialwall.annot');
rh_targ_mesh = CBIG_read_fslr_surface('rh', 'fs_LR_32k','inflated','medialwall.annot');

pos = [lh_targ_mesh.vertices rh_targ_mesh.vertices];
nei_list = lh_targ_mesh.vertexNbors;
nei_s = sort(nei_list, 1, 'ascend');
nei_list = rh_targ_mesh.vertexNbors;
nei_list(nei_list ~= 0) = nei_list(nei_list~= 0) + 32492;
nei_r = sort(nei_list, 1, 'ascend');
nei_v = [nei_s nei_r];

new_label_t = zeros(64984, 1);

x_r_direc = [-0.3861;
    0.7778;
   -1.2339];


for parceli = 1:size(network_match, 1)
    index = find(group_label == parceli);
    if isempty(index) || isempty(network_match{parceli})
        continue;
    end

    if max(size(network_match{parceli})) == 0
        continue;
    elseif max(size(network_match{parceli})) == 1

        new_label_t(index) = network_match{parceli}(1);
        continue;
    end

    

    net_l = network_match{parceli};
    
    v_last_idx = index;
    
    flg_v_cnt = round(length(index) / length(net_l));

    clear x_r_direc
    
    for  ni = 1:length(net_l)-1
        flg_b = 0;
        visit_list = [];
        gflg = 0;
        for vi = 1:length(v_last_idx)
            if length(intersect(nei_v(:, vi), index)) < 6 && isempty(intersect(nei_v(:, vi), visit_list))
                k = vi;gflg = 1;break;
            end
        end
        if gflg == 0
            flg_b = 1;
        end
        
        v_tmp = v_last_idx(k);
        v_last_idx(k) = [];
        new_label_t(v_tmp) = net_l(ni);
        visit_list(end+1) = v_tmp;
    
        while ~flg_b
    
            neibor_tmp = nei_v(:, v_tmp);
            neibor_tmp = intersect(neibor_tmp, index);
    
            
            if ~exist('x_r_direc', 'var')
                for kki = 1:length(neibor_tmp)
                    v_next_tmp = neibor_tmp(kki);
                    if length(intersect(nei_v(neibor_tmp(kki)), index)) == 6
                        break;
                    end
                end
                x_r_direc = pos(:, v_next_tmp) - pos(:, v_tmp);
            else
                v_next_tmp = v_tmp;
            end

            
            new_label_t(v_next_tmp) = net_l(ni);
            v_last_idx = setdiff(v_last_idx,  v_next_tmp);
            visit_list(end+1) =   v_next_tmp;
            visit_list = unique(visit_list);
            
            log_visit = 0;
            while true
                if length(visit_list) > flg_v_cnt || log_visit == length(visit_list)
                    flg_b = 1;
                    break;
                end
                log_visit = length(visit_list);
                neibor_tmp = nei_v(:, v_next_tmp);
                neibor_tmp = intersect(neibor_tmp, index);
      
                x_r_direc_tmp = pos(:, neibor_tmp) - pos(:, v_next_tmp);
                x_r_direc_var = sum(abs(x_r_direc_tmp - x_r_direc), 1);
    
                [ia, ib] = min(x_r_direc_var);
    
                
                if length(intersect(nei_v(:, neibor_tmp(ib)), visit_list)) > 1 || length(intersect(nei_v(:, neibor_tmp(ib)), index)) < 6 
                    break;
                end
                v_next_tmp = neibor_tmp(ib);
                if isempty(intersect(v_last_idx, v_next_tmp))
                    break;
                end
                
                new_label_t(v_next_tmp) = net_l(ni);
                v_last_idx = setdiff(v_last_idx,  v_next_tmp);
                visit_list(end+1) = v_next_tmp;
                visit_list = unique(visit_list);
                
            end
    
            
            flg_t = 0;
            for vi = 1:length(v_last_idx)
                if length(intersect(nei_v(:, v_last_idx(vi)), index)) < 6 && isempty(intersect(nei_v(:, v_last_idx(vi)), visit_list)) 
                    k = vi;
                    flg_t = 1;
                    break;
                end
            end
            if flg_t == 0
                break;
            end
            v_tmp = v_last_idx(k);
        
        end
    end
    new_label_t(v_last_idx) = net_l(end);
end