clc;clear;

addpath('/home/jinlong/VDisk1/Jinlong/2024_homologous_parcellation/1code/HIP2024/example/');
HFP_config_prepare

format_datapath_alls6 = [ getenv('DATASET_DIR') '/HCD/Package_1188937_hcpdrestingRecommand/fmriprep/postproc/91k_AP_original_1.6/data/%s/rfMRI_REST_all_Atlas_MSMAll_hp0_clean_regress_smooth6_all.dtseries.nii'];
format_datapath_sess = [ getenv('DATASET_DIR') '/HCD/Package_1188937_hcpdrestingRecommand/fmriprep/postproc/91k_AP_original_1.6/data/%s/rfMRI_REST_all_Atlas_MSMAll_hp0_clean_regress_smooth6_all_sess%d.dtseries.nii'];

load([getenv('OUT_RET_DIR') '/IHFP_ret/homologous_connectional_homogeneity.mat'], 'new_sub_list', 'new_age_list');

load([getenv('OUT_RET_DIR') '/IHFP_ret/IHFP_group_network_large_orig_homo.mat'], 'age_group_label_gwMRF', 'age_independent_group_label');

load([  '/home/jinlong/VDisk1/Jinlong/IHFP/resources//hcd_age_group/age_group_split_index.mat']);

down_cifti_s = ft_read_cifti('/home/jinlong/VDisk2/Jinlong/resources/fslr_downsample_900mesh_parcellation.dlabel.nii');



age_g_flag = zeros(22, 6);
age_g_flag([8 9], 1) = 1;
age_g_flag([10 11 12], 2) = 1;
age_g_flag([13 14 15], 3) = 1;
age_g_flag([16 17 18], 4) = 1;
age_g_flag([19 20 21], 5) = 1;
age_g_flag(:, 6) = 1;


age_group_corr = cell(6, 2);

age_group_l = cell(5, 1);
for ai = 1:5
    tmp = unique(age_group_label_gwMRF(:, ai));
    tmp(tmp < 1) = [];
    age_group_l{ai} = tmp;
end

age_independent_l = unique(age_independent_group_label);
age_independent_l(age_independent_l < 1)  = [];

age_group_l{6} = age_independent_l;

age_group_all_label = age_group_label_gwMRF;
age_group_all_label(:, 6) = age_independent_group_label;

for sesi = 1:2
    for subi = 1:length(new_sub_list)
            f_l = sprintf(format_datapath_sess, new_sub_list{subi}, sesi);
            f_d = ft_read_cifti(f_l);
            f_t = f_d.dtseries;
            for ai = 1:6
                if age_g_flag(floor(new_age_list(subi)), ai) == 0
                    continue;
                end
                disp(['readdata sub: ' num2str(subi) '/' num2str(length(new_sub_list))]);
                group_1 = age_group_all_label(:, ai);
                group_list = age_group_l{ai};
                new_dt = zeros(length(group_list), size(f_t, 2));
                for ki = 1:length(group_list)
                    idx = find(group_list(ki) == group_1);
                    new_dt(ki, :) = nanmean(f_t(idx, :), 1);
                end
                if min(size(age_group_corr{ai, sesi})) == 0
                    age_group_corr{ai, sesi} = {new_dt};
                else
                    age_group_corr{ai, sesi}{end+1} = new_dt;
                end
            end
            
    end
end

age_group_fc_m = cell(6, 2);
for sesi = 1:2
    for ai = 1:6
        tmp1 = {};
        for ki = 1:length(age_group_corr{ai, sesi})
            disp(['compute roi-wise rsfc group: ' num2str(ai) ' sub: ' num2str(ki) '/' num2str(length(age_group_corr{ai}))]);
            tt = paircorr_mod(age_group_corr{ai, sesi}{ki}');
            tmp1{ki} = tt;
        end
        age_group_fc_m{ai, sesi} = tmp1;
    end
end


age_group_mean_corr_intra = {};

for ai = 1:6
    tmp = zeros(length(age_group_l{ai, 1}), length(age_group_corr{ai, 1}));
    nc1 = 0;
    for si = 1:length(age_group_corr{ai, 1})
        disp(['intra conn  group: ' num2str(ai) '  sub: ' num2str(si) '/' num2str(length(age_group_corr{ai, 1}))]);
        
            corr_tmp = paircorr_mod(age_group_fc_m{ai, 1}{si}', age_group_fc_m{ai, 2}{si}');
            nc1 = nc1 + 1;
            tmp(:, nc1) = diag(corr_tmp);
       
    end
    age_group_mean_corr_intra{ai} = tmp;
end


age_group_mean_corr_inter = cell(6,2);
inter_corr_idx = zeros(length(age_group_corr{ai, 1}));
for ai = 1:6
    tmp = zeros(length(age_group_l{ai, 1}), (length(age_group_corr{ai, 1})-1)*length(age_group_corr{ai, 1})/2);
    nc1 = 0;
    for si = 1:length(age_group_corr{ai, 1})-1
        
        for sj = si+1:length(age_group_corr{ai, 1})
            disp(['inter conn  group: ' num2str(ai) '  sub: ' num2str(si) ' - ' num2str(sj) '  /' num2str(length(age_group_corr{ai, 1}))]);
            corr_tmp = paircorr_mod(age_group_fc_m{ai, 1}{si}', age_group_fc_m{ai, 1}{sj}');
            nc1 = nc1 + 1;
            tmp(:, nc1) = diag(corr_tmp);
            inter_corr_idx(si, sj) = nc1;
            inter_corr_idx(sj, si) = nc1;
        end
    end
    age_group_mean_corr_inter{ai, 1} = tmp;
end

inter_corr_idx = zeros(length(age_group_corr{ai, 1}));
for ai = 1:6
    tmp = zeros(length(age_group_l{ai, 1}), (length(age_group_corr{ai, 1})-1)*length(age_group_corr{ai, 1})/2);
    nc1 = 0;
    for si = 1:length(age_group_corr{ai, 1})-1
        
        for sj = si+1:length(age_group_corr{ai, 1})
            disp(['inter conn  group: ' num2str(ai) '  sub: ' num2str(si) ' - ' num2str(sj) '  /' num2str(length(age_group_corr{ai, 1}))]);
            corr_tmp = paircorr_mod(age_group_fc_m{ai, 2}{si}', age_group_fc_m{ai, 2}{sj}');
            nc1 = nc1 + 1;
            tmp(:, nc1) = diag(corr_tmp);
%             inter_corr_idx(si, sj) = nc1;
%             inter_corr_idx(sj, si) = nc1;
        end
    end
    age_group_mean_corr_inter{ai, 2} = tmp;
end





% Invert the similarity matrix R (1 - R_i)
inter_variability_list = zeros(length(age_group_l{1}), 6);
new_label = zeros(96854, 6);
for ai = 1:6

    R_invert = 1 - mean(age_group_mean_corr_inter{ai}, 2);

    N = 1 - mean(age_group_mean_corr_intra{ai}, 2);
    Y = R_invert;
    X = N;  


    X_mat = [ones(length(X), 1), [1:length(X)]']; 
    b = X_mat \ Y;  

    V = (Y - X_mat * b);  
    inter_variability_list(:, ai) = V;

    gr1 = age_group_all_label(:, ai);
    bound_label = draw_label_boundary(gr1, 'fs_LR_32k');
    gr_list = age_group_l{ai};
    for pppi = 1:length(gr_list)
        idx = find(bound_label == gr_list(pppi));
        if isempty(idx)
            continue;
        end
        new_label(idx, ai) = V(pppi);
    end
end








a = ft_read_cifti('/home/jinlong/VDisk2/Jinlong/resources/test.dtseries.dtseries.nii');
a.dtseries = new_label;
a.dtseries(20248, :) = 0.2;
a.dtseries(20200, :) = -0.14;
a.time = 1:6;
ft_write_cifti([getenv('OUT_RET_DIR') '/intrasubject_connectional_variability'], a, 'parameter', 'dtseries');










