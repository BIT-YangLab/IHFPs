clc;clear;


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

%% compute permutation
prm_metrics = {'corr','COD','predictive_COD','MAE','MAE_norm','MSE','MSE_norm'};


load('/home/jinlong/VDisk2/Jinlong/2024_homologous_parcellation_result/HFIP_result/HFIP_ret/individual_label_list_17_orig_homo.mat', 'sublist_1');

% network-mean feature
fc_mean_network_cell = {[], [], []};
load([pred_src_dir '/' pred_m_dir '.mat']);
% for pi = 1:length(pro_list)
%     eval(['pro_l = ' pro_list{pi} ';']);
%     for subi = 1:length(pred_sub_list)
%         new_fc = reshape(pro_l(subi, :), parcel_n_list(pi), []);
%         [ia, ib] = ismember(pred_sub_list{subi}, sublist_1);
%         if pi == 3
%             sub_assign = all_sub_assign_prob_filter{ib};
%         else
%             eval(['sub_assign = all_sub_assign_' pro_name_list{pi} ';']);
%         end
%         new_fc = new_fc * sub_assign;
%         new_fc = new_fc' * sub_assign;
%         fc_mean_network_cell{pi}(subi, :) = reshape(new_fc, 1, 17*17);
%     end
% end


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
    % FC_network_mean = fc_mean_network_cell{ai}';
    for ti = 1:length(target_dir_list)

        target_v = pred_target(:, ti);
        nan_idx = union(find(isnan(target_v)), find(target_v == 999));
        pred_new_sublist = pred_sub_list;
        pred_new_sublist(nan_idx) = [];
%         FC_network_mean_new = FC_network_mean;
%         FC_network_mean_new(:, nan_idx) = [];
        PFM_score_perm_mean = [];
        for sdi = 1:seed_n
            

            dst_fsm_file_dir = [pred_dst_dir '/' file_dir1 '/' type_dir{ai} '/' target_dir_list{ti} '/seed_' num2str(sdi) '/FSM/FSM_corr.mat'];
            load(dst_fsm_file_dir);

            dst_ret_file_dir = [pred_dst_dir '/' file_dir1 '/' type_dir{ai} '/' target_dir_list{ti} '/seed_' num2str(sdi) '/final_result_example_targets.mat'];
            load(dst_ret_file_dir)
            for i = 1:length(prm_metrics)
                stats_perm.(prm_metrics{i}) = zeros(N_fold,perm_num);
            end

            for fdi = 1:fold_n
                disp(['fold ' num2str(fdi)]);
            end

            y_regress = cell(fold_n,1);
            for fdi = 1:fold_n
                dst_file_dir = [pred_dst_dir '/' file_dir1 '/' type_dir{ai} '/' target_dir_list{ti} '/seed_' num2str(sdi) '/y/fold_' num2str(fdi) '/y_regress_example_targets.mat'];
                load(dst_file_dir, 'y_resid');

                %% load FSM
                

                y_regress{fdi} = y_resid;

                
                opt_lambda = optimal_lambda(fdi,1);
                test_ind = sub_fold(fdi).fold_index;
                train_ind = ~test_ind;
                
                FSM_train = FSM(train_ind,train_ind,:);
                FSM_pred = FSM(test_ind,train_ind,:);
                
                % number of subjects
                num_train_subj = size(FSM_train,1);
                num_test_subj = size(FSM_pred,1);
                
                % Perform training and prediction
                K_r = FSM_train + opt_lambda*eye(num_train_subj);
                % training
                X = ones(num_train_subj,1);
                inv_K_r = inv(K_r);
                beta_pre = (X' * (inv_K_r * X)) \ X' * inv_K_r;
            



                PFM_score_perm = [];
                for pei = 1:perm_num
                    rng(pei+7-1);
                    perm_ind = randperm(length(pred_new_sublist));    
                    y_perm = y_resid(perm_ind);
                    % compute activation
                    y_train_resid = y_perm(train_ind);
                    y_test_resid = y_perm(test_ind);
                    % predict
                    beta = beta_pre * y_train_resid;
                    alpha = inv_K_r * (y_train_resid - X * beta);
                    
                    y_predicted = FSM_pred * alpha + ones(num_test_subj,1) .* beta;
                    for k = 1:length(prm_metrics)
                        stats_perm.(prm_metrics{k})(fdi,pei) = ...
                            CBIG_compute_prediction_acc_and_loss(y_predicted,y_test_resid,prm_metrics{k},y_train_resid);
                    end
                end
                
            
            end
        end

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

function [pred_stats,loss] = CBIG_compute_prediction_acc_and_loss(y_pred,y_test,metric,y_train)

    % [pred_stats,loss] = CBIG_compute_prediction_acc_and_loss(y_pred,y_test,metric,y_train)
    %
    % This function computes the prediction accuracy and loss defined by the
    % metric passed in.
    %
    % Inputs:
    %   - y_pred
    %     Vector of length # test subjects. Prediction of the traget variable of
    %     the test subjects
    %
    %   - y_test
    %     Vector of length # test subjects. Groudtruth of the target variable of
    %     the test subjects
    %
    %   - metric
    %     A string indicating the prediction stats to be computed.
    %     Choose from:
    %       'corr'              - Pearson's correlation;
    %       'COD'               - Coefficient of determination. Defined as
    %                             1-||y_pred-y_test||^2/||mean(y_test)-y_test||^2,
    %                             where y_pred is the prediction of the test data, 
    %                             y_test is the groud truth of the test data, 
    %                             and mean(y_test) is the mean of test data
    %       'predictive_COD'    - Predictive coefficient of determination. Defined as
    %                             1-||y_pred-y_test||^2/||mean(y_train)-y_test||^2,
    %                             where y_pred is the prediction of the test data, 
    %                             y_test is the groud truth of the test data, 
    %                             and mean(y_train) is the mean of training data
    %       'MAE'               - mean absolute error
    %       'MAE_norm'          - mean absolute error divided by the standard
    %                             derivation of the target variable of the training set
    %       'MSE'               - mean squared error
    %       'MSE_norm'          - mean squared error divided by the variance
    %                             of the target variable of the training set
    %
    %   - y_train
    %     Vector of length # test subjects. The target variable of the traning 
    %     subjects. Only useful when compute the predictive_COD, MAE_norm, or MSE_norm
    
    % Outputs:
    %   - pred_stats
    %     A scalar. Prediction statistics of the given metric
    %
    %   - loss
    %     A scalar. Prediction loss of the given the metric. The loss is the same 
    %     as pred_stats if the given metric is the smaller the better (e.g. 'MAE','MSE').
    %     The loss is the addictive inverse of the pred_stats if the given metric is
    %     the bigger the better (e.g. 'corr', 'COD').
    %
    % Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
    
    switch metric
        case 'corr'
            pred_stats = CBIG_corr(y_test,y_pred);
            loss = -pred_stats;
        case 'COD'
            ss_res = sum((y_pred - y_test).^2);
            ss_total = sum((y_test - mean(y_test)).^2);
            pred_stats = 1-ss_res/ss_total;
            loss = -pred_stats;
        case 'predictive_COD'
            ss_res = sum((y_pred - y_test).^2);
            ss_total = sum((y_test - mean(y_train)).^2);
            pred_stats = 1-ss_res/ss_total;
            loss = -pred_stats;
        case 'MAE'
            pred_stats = mean(abs(y_pred-y_test));
            loss = pred_stats;
        case 'MAE_norm'
            pred_stats = mean(abs(y_pred-y_test))/std(y_train);
            loss = pred_stats;
        case 'MSE'
            pred_stats = mean((y_pred - y_test).^2);
            loss = pred_stats;
        case 'MSE_norm'
            pred_stats = mean((y_pred - y_test).^2)/var(y_train);
            loss = pred_stats;
        otherwise
            disp(metric);
            error('Unexpected metric');
    end
end







