clc;
clear;


% pfm_perm_cell



pred_dst_dir = '/home/jinlong/VDisk2/Jinlong/2024_homologous_parcellation_result/HFIP_result/HFIP_ret/prediction_ret';
file_dir1 = 'homologous_rsfc_new_2';
type_dir = {'pred_rsfc_large', 'pred_rsfc_orig', 'pred_rsfc_homo'};

pred_src_dir = '/home/jinlong/VDisk2/Jinlong/2024_homologous_parcellation_result/HFIP_result/prediction_ret';
pred_m_dir = 'prediction_rsfc';

load([pred_src_dir '/prediction_homogeneity.mat'], 'pred_klist', 'pred_sub_list');
target_dir_list = pred_klist;


seed_n = 20;
fold_n = 10;
perm_num = 2000;
load('/home/jinlong/VDisk2/Jinlong/2024_homologous_parcellation_result/HFIP_result/HFIP_ret/group_network_large_orig_homo.mat');
pro_list = {'pred_rsfc_large', 'pred_rsfc_orig', 'pred_rsfc_homo'};
pro_name_list = {'large', 'orig', 'homo'};
parcel_n_list = [34, 400, 400];

load('/home/jinlong/VDisk2/Jinlong/2024_homologous_parcellation_result/HFIP_result/HFIP_ret/individual_label_list_17_orig_homo.mat', 'sublist_1');

% network-mean feature
fc_mean_network_cell = {[], [], []};



load([pred_dst_dir '/' file_dir1 '/no_relative_10_fold_sub_list.mat']);
N_fold = length(sub_fold);
folds = cell(N_fold,1);
for i = 1:N_fold
    test_ind = sub_fold(i).fold_index;
    train_ind = ~test_ind;
    folds{i} = train_ind;
end




pfm_perm_cell = cell(3, 8);
for ai = 1:3
    FC_network_mean = fc_mean_network_cell{ai}';
    for ti = 1:length(target_dir_list)

        target_v = pred_target(:, ti);
        nan_idx = union(find(isnan(target_v)), find(target_v == 999));
        pred_new_sublist = pred_sub_list;
        pred_new_sublist(nan_idx) = [];
        FC_network_mean_new = FC_network_mean;
        FC_network_mean_new(:, nan_idx) = [];
        PFM_score_perm_mean = [];
        for sdi = 1:seed_n
            y_regress = cell(fold_n,1);
            for fdi = 1:fold_n
                dst_file_dir = [pred_dst_dir '/' file_dir1 '/' type_dir{ai} '/' target_dir_list{ti} '/seed_' num2str(sdi) '/y/fold_' num2str(fdi) '/y_regress_example_targets.mat'];
                load(dst_file_dir, 'y_resid');
                y_regress{fdi} = y_resid;
            end
            PFM_score_perm = [];
            for pei = 1:perm_num
                rng(pei+7-1);
                perm_ind = randperm(length(pred_new_sublist));    
                % compute activation
                activation_score_fold = zeros(17*17,1,N_fold);
                for fdj = 1:fold_n
                    train_ind = folds{fdj};
                    y_resid = y_regress{fdj};
                    y_resid = y_resid(perm_ind,:);
                    y_resid = y_resid(train_ind,:);
                    features_train = FC_network_mean_new(:,train_ind);
                    curr_cov = CBIG_TRBPC_compute_covariance(features_train',y_resid);
                    activation_score_fold(:,:,fdj) = bsxfun(@times,curr_cov,1./std(y_resid));
                end
                
                PFM_score_perm(:,pei) = mean(activation_score_fold,3);
            end
            if sdi == 1
                PFM_score_perm_mean =  PFM_score_perm;
            else
                PFM_score_perm_mean =  PFM_score_perm_mean + PFM_score_perm;
            end
        end
        PFM_score_perm_mean = PFM_score_perm_mean ./ seed_n;
        pfm_perm_cell{ai, ti} = PFM_score_perm_mean;  
    end

end






function cov_xy = CBIG_TRBPC_compute_covariance(x,y)

% cov_xy = CBIG_TRBPC_compute_covariance(x,y)
%
% x: N*K1 matrix, N is # observation, K is # variables
% y: N*K2 matrix, N is # observation, K is # variables
% cov_xy: K1*K2 covariance matrix
%
% Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cov_xy = bsxfun(@minus,x,mean(x))'*bsxfun(@minus,y,mean(y))/(size(x,1)-1);

end












