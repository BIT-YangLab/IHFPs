


load([ getenv('OUT_RET_DIR') '/IHFP_ret/individual_label_list_17_orig_homo.mat']);

a = load([ getenv('OUT_RET_DIR') '/IHFP_ret/IHFP_sess1.mat']);

sublist_all = sublist_1;


load([getenv('RESOURCE_DIR') '/fs_LR_32k_no_medial_wall_index.mat']);  %no_medial_wall_index
ref_group = load([getenv('RESOURCE_DIR')  '/HCP_40sub_1000iter_17nets_cen_sm4.mat']);
label_ref = [ref_group.lh_labels; ref_group.rh_labels];




format_data_path = '/home/jinlong/VDisk2/Jinlong/data//HCD/Package_1188937_hcpdrestingRecommand/fmriprep/postproc/91k_AP_original_1.6/data/%s/rfMRI_REST_all_Atlas_MSMAll_hp0_clean_regress_smooth6_all.dtseries.nii';
network_n = 400;

segregation_vlist = zeros(length(sublist_all), network_n);
segregation_parcel_list = cell(network_n, 1);
large_scale_segregation_list = zeros(length(sublist_all), network_n);




mean_fc_list = zeros(length(sublist_all), 1);
mean_std_list = zeros(length(sublist_all), 1);
fc_parcel_list_within = zeros(size(sublist_all, 2), network_n);
std_parcel_list_within = zeros(size(sublist_all, 2), network_n);

fc_parcel_list_between = zeros(size(sublist_all, 2), network_n);
std_parcel_list_between = zeros(size(sublist_all, 2), network_n);

fc_parcel_list_all = zeros(size(sublist_all, 2), network_n);
std_parcel_list_all = zeros(size(sublist_all, 2), network_n);


for subi = 1:length(sublist_all)
    subname = sublist_all{subi}; tic;
    disp(['process subi:' num2str(subi) '/' num2str(length(sublist_all)) ' - ' subname])
    cifti_f = sprintf(format_data_path, subname);
    if ~exist(cifti_f)
        continue;
    end

   

    label_ind = label_list{3}(:, subi);
    label_ind(isnan(label_ind)) = 0;
    
    [homo_ci, std_ci ] = compute_parcellation_mean_fc(label_ind, {cifti_f});


    label_t_list = unique(label_ind);
    label_t_list(label_t_list < 1) = [];
    label_t_list(isnan(label_t_list)) = [];
    fc_parcel_list_within(subi, label_t_list) = homo_ci(:, 1);
    std_parcel_list_within(subi, label_t_list) = std_ci(:, 1);

    fc_parcel_list_between(subi, label_t_list) = homo_ci(:, 2);
    std_parcel_list_between(subi, label_t_list) = std_ci(:, 2);

    fc_parcel_list_all(subi, label_t_list) = homo_ci(:, 3);
    std_parcel_list_all(subi, label_t_list) = std_ci(:, 3);

    toc;
end

save([getenv('OUT_RET_DIR')  '/homologous_global_mean_fc.mat'], 'fc_parcel_list_within', 'fc_parcel_list_between', 'fc_parcel_list_all');


