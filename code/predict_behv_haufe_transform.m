clc;clear;

src_dir = '/home/jinlong/VDisk2/Jinlong/2024_homologous_parcellation_result/HFIP_result/prediction_ret/prediction_homogeneity.mat';
ret_dir = '/home/jinlong/VDisk2/Jinlong/2024_homologous_parcellation_result/HFIP_result/HFIP_ret/prediction_ret/homologous_homogeneity_sa/pred_homogeneity_homo/seed_%d/final_result_example_targets.mat';
fold_dir = '/home/jinlong/VDisk2/Jinlong/2024_homologous_parcellation_result/HFIP_result/HFIP_ret/prediction_ret/homologous_homogeneity_sa/no_relative_20_fold_sub_list.mat';

load(src_dir);
load(fold_dir);



PFM_all = cell(2, 1);
PFM_all_folds = zeros(400, 20);

for te = 1:2
    for fi = 1:20
        load(sprintf(ret_dir, fi));
        features_train = pred_homogeneity_homo;
        fold_idx = sub_fold(fi).fold_index;
        train_idx = fold_idx == 0;
        test_idx = fold_idx == 1;
        y_train_fe = y_predict_concat(train_idx, te);
        y_train_pred = y_pred_train{fi}(:, te);
        nan_idx = union(find(isnan(y_train_fe)), find(isnan(y_train_pred)));
        y_train_fe(nan_idx) = [];
        y_train_pred(nan_idx) = [];
        features_train = features_train(train_idx, :);
        features_train(nan_idx, :) = [];
        PFM_all_folds(:,fi) = CBIG_PFM_cov_matrix(features_train, y_train_pred)/std(y_train_pred);
    
    end
    PFM_all{te} = PFM_all_folds;


end




