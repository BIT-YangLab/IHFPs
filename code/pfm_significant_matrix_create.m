

pred_dst_dir = [ getenv('OUT_RET_DIR') '/HFIP_result/prediction_ret'];
file_dir1 = 'pred_rsfc_between';
type_dir = {'pred_rsfc_orig', 'pred_rsfc_homo'};

if ~exist([pred_dst_dir '/permutation_ret'])
    mkdir([pred_dst_dir '/permutation_ret']);
end

seed_n = 20;
fold_n = 10;
perm_num = 200;
perm_nc = 500;
load([ getenv('OUT_RET_DIR') '/HFIP_result/HFIP_ret/group_network_large_orig_homo.mat']);
pro_list = {'pred_rsfc_orig', 'pred_rsfc_homo'};
pro_name_list = {'orig', 'homo'};
parcel_n_list = [400, 400];

load('/HFIP_result/HFIP_ret/individual_label_list_17_orig_homo.mat', 'sublist_1');

% network-mean feature
fc_mean_network_cell = {[], [], []};
load( [ getenv('OUT_RET_DIR') '/HFIP_result/prediction_ret/pred_rsfc_between.mat']);
target_dir_list = pred_klist;


for pi = 1:length(pro_list)
    eval(['pro_l = ' pro_list{pi} ';']);
    eval(['resample_idx = resample_idx_' pro_name_list{pi} ';']);
    for subi = 1:length(pred_sub_list)
        new_fc = zeros(400, 400);
        new_fc(resample_idx) = pro_l(subi, :);
        
        [ia, ib] = ismember(pred_sub_list{subi}, sublist_1);
        if pi == 2
            sub_assign = all_sub_assign_prob_filter{ib};
        else
            eval(['sub_assign = all_sub_assign_' pro_name_list{pi} ';']);
        end
        new_fc = new_fc * sub_assign;
        new_fc = new_fc' * sub_assign;
        new_fc1 = new_fc(eye(size(new_fc, 1)) == 0);
        fc_mean_network_cell{pi}(subi, :) = reshape(new_fc1, 1, []);
    end
end


load([pred_dst_dir '/' file_dir1 '/no_relative_10_fold_sub_list.mat']);
N_fold = length(sub_fold);
folds = cell(N_fold,1);
for i = 1:N_fold
    test_ind = sub_fold(i).fold_index;
    train_ind = ~test_ind;
    folds{i} = train_ind;
end




pfm_perm_cell = cell(2, 8);
for ai = 1:2
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
            for pei = 1:1

  
                % compute activation
                activation_score_fold = zeros(size(FC_network_mean, 1),1,N_fold);
                for fdj = 1:fold_n
                    train_ind = folds{fdj};
                    if ti == 4
                        train_ind = train_ind(1:end-1);
                    end
                    y_resid = y_regress{fdj};

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
save([ pred_dst_dir '/permutation_ret/permutation_pfm_homo_raw.mat'], 'pfm_perm_cell');
toc;




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







