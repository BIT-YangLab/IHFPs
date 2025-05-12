
function group_fine2network_profile(out_dir, group_homo_atlas)
% DESCRIPTION:
%   Calculate similarity profile between parcel and network based on group atlas
%
% USAGE: 
%   group_fine2network_profile(out_dir, group_homo_atlas)
%
% Inputs:   out_dir: directory of result
%
%           group_homo_atlas: group large-scale atlas
%
% Outputs:  
% Written by Jinlong Li, Guoyuan Yang
    seed_mesh = 'fs_LR_900';
    targ_mesh = 'fs_LR_32k';

    network_n = 17;
    num_sess = 1;

    load([getenv("RESOURCE_DIR") '/subject_info.mat']);

    sub_list_all_file = [ getenv('RESOURCE_DIR') '/hcd_age_group/hcd_sub_qc_new.txt'];
    sublist_all = sub_info.src_subject_id;
    writecell(sublist_all, sub_list_all_file);

    load([getenv('RESOURCE_DIR') '/fs_LR_32k_no_medial_wall_index.mat']);  %no_medial_wall_index
    if network_n == 7
        ref_group = load([getenv('RESOURCE_DIR')  '/HNU_example_clusters07_scrub_cifti_post_smooth_standard.mat']);
    else
        ref_group = load([getenv('RESOURCE_DIR')  '/HCP_40sub_1000iter_17nets_cen_sm4.mat']);
    end
    label_ref = [ref_group.lh_labels; ref_group.rh_labels];

    
    mkdir([out_dir '/all_profile']);
    temp_txt_profile = [out_dir '/test_profile.txt'];




    for subi = 1:length(sublist_all)
        subname = sublist_all{subi};
        for sess = 1:num_sess
            cifti_file = sub_info.file_all{subi};
            fid = fopen(temp_txt_profile, 'w');
            fprintf(fid, '%s\n', cifti_file);
            fclose(fid);
            profile_file = [ out_dir '/all_profile/' subname '_sess' num2str(sess) '.mat' ];
            profile_file_raw = [ out_dir '/all_profile/' subname '_sess' num2str(sess) '_raw.mat' ];
            CBIG_ComputeCorrelationProfile(seed_mesh,targ_mesh, profile_file, 'NONE', '0.1', temp_txt_profile, 'NONE', 'NONE', '0', profile_file_raw);
        end

    end


    num_data = 0;
    for subki = 1:length(sublist_all)

        subname = (sublist_all{subki});
        profile_file = [out_dir '/all_profile/' subname '_sess' num2str(1) '.mat' ];
        if exist(profile_file)
            num_data = num_data + 1;
            
            profile_data = load(profile_file);

            if(num_data == 1)
                avg_profile = profile_data;
            else
                avg_profile.profile_mat = avg_profile.profile_mat + profile_data.profile_mat;
            end
        end
    end


    group_cnt = 0;
    large_overlap_cell = cell(num_sess,1);
    group_cnt = group_cnt + 1;
    profile_mat = avg_profile.profile_mat./num_data;
    save([out_dir '/avg_profile1.mat'],'profile_mat','-v7.3');
    CBIG_VonmisesSeriesClustering_fix_bessel_randnum_bsxfun(targ_mesh, 'NONE', network_n, [out_dir '/group_large1.mat'], ...
    [out_dir '/avg_profile1.mat'], 'NONE', 0, 1000, 0, 100, 1);

    a1 = load([out_dir '/group_large1.mat']);
    label = [a1.lh_labels; a1.rh_labels];

    [output, assign, cost, dice_overlap] = CBIG_HungarianClusterMatch(label_ref, label);
    large_label_ref = output;

    large_voxel_cell = output;


    for sess = 1:num_sess
        
        labels = group_homo_atlas;

        label_l = unique(labels);
        label_l(label_l<1) = [];

        group_fine_profile_parcel = zeros(size(label_l, 1), size(profile_mat, 2));
        
        new_profile_voxel = zeros(size(profile_mat));
        for i = 1:length(label_l)
            index = find(labels == label_l(i));
            
            group_fine_profile_parcel(i, :) = nanmean(profile_mat(index, :), 1);
        end
        
        label_list = unique(labels);
        label_list(label_list < 1) = [];
        temp_overlap_label = zeros(64984, 1);
        temp_network_match = zeros(length(label_list), 1);
        for li = 1:length(label_list)
            index = find(labels == label_list(li));
            n_t_list = large_label_ref(index);
            n_t_list(n_t_list < 1) = [];
            n_mode = mode(n_t_list);
            temp_overlap_label(index) = n_mode;
            temp_network_match(li) = n_mode;
        end
        bound_label = draw_label_boundary(labels, 'fs_LR_32k');
        temp_overlap_label(bound_label == 0 ) = 0;
        large_overlap_cell{sess} = temp_overlap_label; 
    end


    group_fine_label = labels;
    group_fine_label_list = label_list;

    group_large_label = large_label_ref;
    group_network_match = temp_network_match;
    group_voxel_profile = avg_profile.profile_mat./num_data;
    grp_network.large_overlap_cell = large_overlap_cell;
    grp_network.large_voxel_cell = large_voxel_cell;

    save([ out_dir '/age_independent_template_profile_17net.mat'], 'group_fine_label', 'group_fine_label_list', ...
        'group_fine_profile_parcel', 'group_large_label', 'group_network_match', 'group_voxel_profile', 'grp_network');
end