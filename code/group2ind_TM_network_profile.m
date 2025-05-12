function group2ind_TM_network_profile(out_dir, group_fine_label, new_label_list_all, all_sub_list)
% DESCRIPTION:
%   Calculate similarity profile between group and individual atlas
%
% USAGE: 
%   group2ind_TM_network_profile(out_dir, group_fine_label, new_label_list_all, all_sub_list)
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

    %% process flg
    flg_have_profile = 1;
    flg_have_raw_profile = 0;


    %% data prepare

    load([getenv("RESOURCE_DIR") '/subject_info.mat']);
    network_n = 17;


    load([getenv('RESOURCE_DIR') '/fs_LR_32k_no_medial_wall_index.mat']);  %no_medial_wall_index
    if network_n == 7
        ref_group = load([getenv('RESOURCE_DIR')  '/HNU_example_clusters07_scrub_cifti_post_smooth_standard.mat']);
    else
        ref_group = load([getenv('RESOURCE_DIR')  '/HCP_40sub_1000iter_17nets_cen_sm4.mat']);
    end



    output_dir = [ out_dir '/network_match7'];
    mkdir([out_dir '/all_profile']);
    temp_txt_profile = [out_dir '/test_profile.txt'];


    % create individual connectome profile and group template profile

    load([ out_dir '/age_independent_template_profile_17net.mat']);


    % each age group label

    downsample_dlabel = ft_read_cifti([ getenv('RESOURCE_DIR') '/fslr_downsample_900mesh_parcellation.dlabel.nii']);

    % individual atlas
    % load([ getenv('OUT_RET_DIR') '/HFIP_ret/homologous_ret_sessall_fill.mat'] );


    group_fine_label_list_L = unique(group_fine_label(1:32492));
    group_fine_label_list_L(group_fine_label_list_L < 1) = [];
    group_fine_label_list_R = unique(group_fine_label(32493:64984));
    group_fine_label_list_R(group_fine_label_list_R < 1) = [];
    group_fine_label_list = unique(group_fine_label); group_fine_label_list(group_fine_label_list < 1) = [];

    group_fine_profile_parcel_L = group_fine_profile_parcel(ismember(group_fine_label_list, group_fine_label_list_L), :);
    group_fine_profile_parcel_R = group_fine_profile_parcel(ismember(group_fine_label_list, group_fine_label_list_R), :);

    new_label_list_c = [];

    for ki = 1:5
        label_temp_list = new_label_list_all{ki};
        new_label_list_c = [new_label_list_c label_temp_list];
    end

    template_matching_ret = cell(length(all_sub_list), 1);
    individual_profile_parcel = cell(length(all_sub_list), 1);

    for subi = 1:length(all_sub_list)
        tic
        subname = all_sub_list{subi};
        disp(['process ' num2str(subi) '/' num2str(length(all_sub_list)) ':  ' subname])

        if flg_have_raw_profile == 0
            cifti_file = sub_info.file_all{subi};
            cifti_s = ft_read_cifti(cifti_file);

            good_TR = ~isnan(cifti_s.dtseries(1, :));

            cifti_data = cifti_s.dtseries(1:64984, logical(good_TR));

            cifti_corrmap = paircorr_mod(cifti_data', cifti_data(downsample_dlabel.dlabel ~= 0, :)');
            
            
            % Remove NaNs (produced if vertices have no data)
            cifti_corrmap(isnan(cifti_corrmap)) = 0;
            % Apply the Fisher tranformation
            cifti_corrmap = single(FisherTransform(cifti_corrmap));
        else
            profile_file_raw = [ out_dir '/all_profile/' subname '_sess' num2str(1) '_raw.mat' ];
            cifti_corrmap = load(profile_file_raw);
            cifti_corrmap = cifti_corrmap.corr_raw_max;
        end

        labels_L = new_label_list_c(1:32492, subi);
        label_list_L = unique(labels_L);
        label_list_L(label_list_L < 1) = [];
        labels_R = new_label_list_c(32493:64984, subi);
        label_list_R = unique(labels_R);
        label_list_R(label_list_R < 1) = [];


        new_profile_parcel_L = zeros(length(label_list_L), size(cifti_corrmap, 2));
        new_profile_parcel_R = zeros(length(label_list_R), size(cifti_corrmap, 2));
        for li1 = 1:length(label_list_L)
            index = find(labels_L == label_list_L(li1));
            new_profile_parcel_L(li1, :) = nanmean(cifti_corrmap(index, :), 1);
        end
        for li1 = 1:length(label_list_R)
            index = find(labels_R == label_list_R(li1));
            new_profile_parcel_R(li1, :) = nanmean(cifti_corrmap(index, :), 1);
        end

        individual_profile_parcel{subi} = [new_profile_parcel_L; new_profile_parcel_R];
        eta_to_template_vox_L = template_matching_eta2(group_fine_profile_parcel_L, new_profile_parcel_L);
        eta_to_template_vox_R = template_matching_eta2(group_fine_profile_parcel_R, new_profile_parcel_R);

        [m_v_L, m_i_L] = maxk(eta_to_template_vox_L, round(size(group_fine_profile_parcel_L, 1)*0.05), 2); 
        [m_v_R, m_i_R] = maxk(eta_to_template_vox_R, round(size(group_fine_profile_parcel_R, 1)*0.05), 2); 
        m_i_R = m_i_R + length(group_fine_label_list_L);
        m_i = [m_i_L; m_i_R];
        network_match = group_network_match(m_i);

        template_matching_ret{subi} = network_match;

        toc
    end

    save([ output_dir '/fine2network_TM_ret_17net.mat'], 'template_matching_ret', 'individual_profile_parcel', '-v7.3');



end





