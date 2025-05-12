


function homologous_ind_match(out_dir, group_fine_label, new_label_list_all, all_sub_list)
% DESCRIPTION:
%   Individualized homologous functional parcellations procedure
%
% USAGE: 
%   homologous_ind_match(out_dir, group_fine_label, new_label_list_all, all_sub_list)
%
% Inputs:   out_dir: directory of result
%
%           group_fine_label: group fine-grained atlas
%           
%           new_label_list_all: IHFPs
%       
%           all_sub_list: list of subject 
%
% Outputs:     
% Written by Jinlong Li, Guoyuan Yang
    HFP_config_prepare
    
    load([ out_dir '/age_independent_template_profile_17net.mat'])


    network_n = 17;

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



    new_label_list_c = [];

    for ki = 1:5
        label_temp_list = new_label_list_all{ki};
        new_label_list_c = [new_label_list_c label_temp_list];
    end


    % template_matching_ret    n_sub*1 cell,   398*20
    load([ out_dir '/fine2network_TM_ret_17net.mat']);

    grp_parcel_n = 400;


    all_sim_assign = zeros(grp_parcel_n, network_n);

    n_sub = length(all_sub_list);
    label_check_flg = zeros(400, 1);

    all_sub_assign_m = zeros(400*network_n, n_sub);

    all_sub_assign_prob = cell(n_sub, 1);

    for subi = 1:n_sub

        labels = new_label_list_c(:, subi);
        label_list = unique(labels);
        label_list(label_list < 1) = [];
        label_check_flg(label_list) = label_check_flg(label_list) + 1;

        match_net_temp = template_matching_ret{subi};

        

        temp_net_assi = zeros(400, network_n);
        temp_net_assi_prob = zeros(400, network_n);

        for neti = 1:network_n
            temp_ns = match_net_temp == neti;
            temp_net_assi_prob(label_list, neti) = sum(temp_ns, 2) ;
        end
        all_sub_assign_prob{subi} = temp_net_assi_prob;


        for tai = 1:size(match_net_temp, 1)
            temp_net_assi_cnt = zeros(1, network_n);
            temp_net_assi_cnt(unique(match_net_temp(tai, :))) = 1;
            all_sim_assign(label_list(tai), :) = all_sim_assign(label_list(tai), :) + temp_net_assi_cnt;
            temp_net_assi(label_list(tai), :) = temp_net_assi_cnt;
        end
        temp_net_assi = reshape(temp_net_assi, [], 1);

        

        all_sub_assign_m(:, subi) = temp_net_assi;
        
    end

    temp_parcel_l = ones(400, network_n);
    temp_network_l = ones(400, network_n);

    temp_network_l = temp_network_l * diag(1:network_n);
    temp_parcel_l = diag(1:400) * temp_parcel_l;
    parcel_network_match_list = [ reshape(temp_parcel_l, [], 1) reshape(temp_network_l, [], 1)];

    
    

    %% age-independent
    temp_all_sub_homo_list = cell(400, 1);
    B3 = sum(all_sub_assign_m, 2) ./ 591;
    B2 = parcel_network_match_list;
    B1 = all_sub_assign_m;
    [sv, si] = sort(B3, 'descend');
    B3_a = B3(si);
    B1_a = B1(si, :);
    B2_a = B2(si, :);
    si = find(B3_a > 0.92);
    B2_b = B2_a(si, :);
    % useless
    si = ismember(B2_a(:, 1), B2_b(:, 1));
    B1_c = B1_a(~si, :);
    B2_c = B2_a(~si, :);
    B3_c = B3_a(~si, :);
    %=
    for i = 1:size(B2_b, 1)
        if min(size(temp_all_sub_homo_list{B2_b(i, 1)})) == 0
            temp_all_sub_homo_list{B2_b(i, 1)} = B2_b(i, 2);
        else
            temp_all_sub_homo_list{B2_b(i, 1)} = [temp_all_sub_homo_list{B2_b(i, 1)}  B2_b(i, 2)];
        end
    end
    bound_label = draw_label_boundary(group_homo_label1, 'fs_LR_32k');
    age_group_all_sub_label =  network_striped(bound_label, temp_all_sub_homo_list);

    homo_match_ret = [B2_a(si, :) B3_a(si, :)];
    homo_match_system = cell(max(B2_b(:, 1)), 1);
    homo_system_list = cell(max(B2_b(:, 2)), 1);
    for ki = 1:max(B2_b(:, 1))
        homo_match_system{ki} = B2_b(find(B2_b(:, 1) == ki), 2);
    end

    for ki = 1:max(B2_b(:, 2))
        homo_system_list{ki} = B2_b(find(B2_b(:, 2) == ki), 1);
    end



    save([out_dir '/homologous_large_network_parcel_match.mat'], 'B1_a', 'B2_a', 'B3_a',  ...
        'age_group_all_sub_label', 'all_sub_assign_m', 'homo_match_ret', 'homo_match_system', ...
        'homo_system_list', 'all_sub_assign_prob');

end







