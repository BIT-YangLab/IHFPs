function compute_trajectories_variability(src_mat, dst_p_mat, dst_n_mat)
% DESCRIPTION:
%   Calculate the variability of normal trajectories at the area- and network- level, 
% 　which is the average MAD of each point on the horizontal axis of two curves. 
% 　Then, use GLM regression to determine the number of partitions in each network.
%
% USAGE: 
%   compute_trajectories_variability('dev_fc_all_homo.mat', 'dev_fc_all_homo', 'dev_fc_all_homo_mean_network_weight')
%
% Inputs:   src_mat : areal-level functional features
%
%           dst_p_mat: directory about result of areal-level normative trajectories 
%        
%           dst_n_mat: directory about result of network-level normative trajectories
%
%
% Outputs:  
%
% Written by Jinlong Li, Guoyuan Yang

    addpath('/home/jinlong/VDisk1/Jinlong/2024_homologous_parcellation/example/');
    HFP_config_prepare
    out_dir = getenv("OUT_RET_DIR");
    load([ out_dir '/curv_info/' src_mat ])
    
    load([ out_dir '/prediction_ret/dev_fc_all_homo/gamlss_fit_allData_rs_global_fc_7.mat'], 'all_cen_x');
    load([ out_dir '/HFIP_ret/homologous_large_network_parcel_match_17_new2.mat']);
    
    load([ out_dir '/HFIP_ret/group_network_large_orig_homo.mat' ]);
    
    label_lg = unique(group_homo_label1);
    label_lg(label_lg < 1) = [];
    label_lg(isnan(label_lg)) = [];
    label_s_l = zeros(size(label_lg));
    for pi = 1:400
        idx = find(group_homo_label1 == pi);
        if isempty(idx)
            continue;
        end
        label_s_l(pi) = length(idx);
    end
    label_s_l = label_s_l ./ sum(label_s_l);
    
    index = intersect(find(all_cen_x(:, 1) >= 8), find((all_cen_x(:, 1) < 22)));
    
    
    dv_cell = cell(17, 1);
    
    cen_y_list = zeros(length(index), 400);
    dev_y_list = zeros(8000, 400);
    cnt_p_k = zeros(400, 1);
    
    match_pp = d_curv.parcel_list;
    
    
    for pi = 1:length(match_pp)
        at = load(['/home/jinlong/VDisk2/Jinlong/2024_homologous_parcellation_result/HFIP_result/prediction_ret/' dst_p_mat '/gamlss_fit_allData_rs_global_fc_' num2str(pi+1) '.mat']);
        cen_y_list(:, match_pp(pi)) = at.all_cen_y(index, 4);
        dev_y_list(:, match_pp(pi)) = dev_y_list(:, match_pp(pi)) + at.all_velocity;
        cnt_p_k(match_pp(pi)) = cnt_p_k(match_pp(pi)) + 1;
    end
    
    
    dv_list = [];
    parcel_n_l = [];
    
    for neti = 1:17
        i2 = homo_system_list{neti};
        if isempty(i2)
            continue;
        end
    
        curv_list = zeros(length(index), length(i2));
        temp_prob =  [];
        for i2_i = 1:length(i2)
            idx1 = intersect(find(homo_match_ret(:, 2) == neti), find(homo_match_ret(:, 1) == i2(i2_i)));
            curv_list(:, i2_i) = cen_y_list(:, i2(i2_i));
            temp_prob(end+1) = homo_match_ret(idx1, 3);
        end
    
        gam_nc = neti + 1;
        if neti > 8
            gam_nc = neti;
        end
    
        bt = load(['/home/jinlong/VDisk2/Jinlong/2024_homologous_parcellation_result/HFIP_result/prediction_ret/' dst_n_mat '/gamlss_fit_allData_rs_global_fc_' num2str(gam_nc) '.mat']);
        curv_n = bt.all_cen_y(index, 4);
        dv = nanmean(abs(curv_n - curv_list), 1).*temp_prob; dv = dv';
        dv_cell{neti} = dv;
    
        dv_list = [dv_list; dv];
        
        parcel_n_l = [parcel_n_l; label_s_l(i2)];
    end
    
    uidx = unique(homo_match_ret(:, 1));
    
    [ia, ib] = ismember(uidx, homo_match_ret(:, 1));
    dv_list1 = dv_list(ib);
    parcel_n_l1 = parcel_n_l(ib);
    
    [r_res, beta, ~, ~] = CBIG_glm_regress_matrix(dv_list1, parcel_n_l1);
    [ia, ib] = ismember(homo_match_ret(:, 1), uidx);
    A_pred = r_res(ib);
    
    start_idx = 1; end_idx = 0;
    for neti = 1:17
        if isempty(dv_cell{neti})
            continue;
        end
        i2 = homo_system_list{neti};
        end_idx = start_idx + length(i2) - 1;
        dv_cell{neti} = A_pred(start_idx:end_idx);
        start_idx = end_idx + 1;
    end
    save([out_dir '/prediction_ret/gamlss_info_ret/diff_cell_' dst_p_mat '.mat'], "dv_cell");
end

    