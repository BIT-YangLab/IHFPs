function predict_behv_model_run(temp_multi_stem, process_queue, feature_mat_path_dir)
% modify based on exe_example_KRR_wrapper.m, written by Ruby Kong, CBIG
%
% This scripts is a wrapper function to run the LRR_fracfridge workflow using example
% data for two resolutions 100 and 200. 
%
% Inputs:
%   - out_dir
%     The output directory where the results will be saved.
%
%  output (saved)
%   - acc_metric_train - training accuracy
%   - acc_corr_test - test accuracies (given in correlation)
%   - y_predict - predicted target values
%   - optimal_statistics - a cell array of size equal to the number of
%   folds. Each element in the cell array is a structure storing the
%   accuracies of each possible accuracy metric (eg. corr, MAE, etc).
% 
% Written by Ruby Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
% CBIG_KRR_test_cv_allparams_LITE bin mod
% CBIG_KRR_innerloop_cv_allparams_LITE bin mod

addpath('/home/jinlong/VDisk1/Jinlong/2024_homologous_parcellation/example/');
HFP_config_prepare
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');


% % temp_multi_stem = 'orig_homogeneity';
out_dir = [ getenv('OUT_RET_DIR') '/HFIP_ret/prediction_ret/' temp_multi_stem]; 

if ~exist(out_dir)
    mkdir(out_dir);
end


num_folds = 10;
% seed = 1;

delimiter = ',';



% % process_queue = {'pred_homogeneity'};
feature_queue = { 'all_fc_list', 'all_fc_list_grp1', 'all_fc_list_grp2', 'all_fc_list_grp3'};
% target_queue = {'nih_fluidcogcomp_unadjusted'};%'er_rtcr', 'er_ang', 'er_fear', 'er_hap', 'er_noe', 'er_sad'};

% sub_fold_path = fullfile(DATA_DIR, 'no_relative_20_fold_sub_list.mat');
% feature_mat_path = load([getenv('OUT_RET_DIR') '/prediction_ret/prediction_homogeneity_orignal_kong2022.mat']);%matfile(fullfile(DATA_DIR, '13_15_age_label_all_prediction.mat'));
feature_mat_path = load(feature_mat_path_dir);%matfile(fullfile(DATA_DIR, '13_15_age_label_all_prediction.mat'));

target_n = length(feature_mat_path.pred_klist);

for ti = 1:target_n

   target_v = feature_mat_path.pred_target(:, ti);
   nan_idx = union(find(isnan(target_v)), find(target_v == 999));
   

   
   sub_list = feature_mat_path.pred_sub_list;
   sub_list(nan_idx) = [];

   subject_new_path = [out_dir '/new_temp_sublist.txt'];
   writecell(sub_list, subject_new_path );
   for seed = 1:20
       sub_fold = CBIG_cross_validation_data_split( subject_new_path, 'none', '', '', num_folds, seed, out_dir, delimiter );
       for pr_i = 1:length(process_queue)
           
               disp(['process ' process_queue{pr_i}]);

               out_file = fullfile(out_dir, process_queue{pr_i}, feature_mat_path.pred_klist{ti} , ['seed_' num2str(seed)], 'final_result_example_targets.mat');
               
               % Check if the output file already exists
               if (~exist(out_file))
                   % Set up parameters
                   % Set up cross-validation folds
   %                 sub_fold_path = matfile(sub_fold_path);
                   params.sub_fold = sub_fold;

                   % Set up feature matrix
                   
                   % feature_names = fieldnames(feature_mat_path);
                   info_struct = feature_mat_path.(process_queue{pr_i});
                   info_struct(nan_idx, :) = [];
                   params.feature_mat = info_struct;
                   params.feature_mat = params.feature_mat(:, :)';
                   temp_fe_mat = fullfile(out_dir, 'feature_temp.mat');
                   temp_fe = params.feature_mat;
                   save(temp_fe_mat, "temp_fe");

                   % Set up covariates matrix
                   % covariates_path = matfile(fullfile(DATA_DIR, 'regressors.mat'));
                   
% %                    params.covariates = feature_mat_path.pred_fd_list;%[feature_mat_path.pred_age_list./12 feature_mat_path.pred_sex_list feature_mat_path.pred_fd_list];%[info_struct.all_sex_list(index_flt) info_struct.all_age_list(index_flt)];
                   params.covariates = [feature_mat_path.pred_age_list./12 feature_mat_path.pred_sex_list feature_mat_path.pred_fd_list];%[info_struct.all_sex_list(index_flt) info_struct.all_age_list(index_flt)];
                   params.covariates(nan_idx, :) = [];

                   %% TODO
                   % Set up target matrix
%                    y_path = matfile(fullfile(DATA_DIR, 'targets.mat'));
                   params.y = feature_mat_path.pred_target(:, ti);
                   params.y(nan_idx, :) = [];
% %                    params.y = [feature_mat_path.pred_age_list./12 feature_mat_path.pred_sex_list];

                   % Set up output directory
                   params.outdir = fullfile(out_dir, process_queue{pr_i}, feature_mat_path.pred_klist{ti}, ['seed_' num2str(seed)]);

                   % Set up output file name
                   params.outstem = 'example_targets';

                   % Set up parameters for LRR_fracridge: number of innerloop cross-validation folds
                   params.num_inner_folds = 10;

                   % Set up metric used for tuning hyperparameter
                   params.metric = 'predictive_COD';

                   % Set up bias flag
                   params.with_bias = 1;

                   % Set up kernel type
                   params.ker_param.type = 'corr';
                   params.ker_param.scale = NaN;

                   % Set up regularization parameter searching range
                   params.lambda_set = [ 0 0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 0.4 0.7...
                       1 1.5 2 2.5 3 3.5 4 5 10 15 20];

                   % Call KRR workflow
                   CBIG_KRR_workflow_LITE(params);
                   
                   for scorei = 1:size(params.y, 2)
                       CBIG_compute_singleKRR_PFM(temp_fe_mat, params.outdir, ...
                           fullfile(out_dir, ['no_relative_' num2str(num_folds) '_fold_sub_list.mat']), ...
                           scorei, [params.outstem '.mat'], params.outdir);
                   end

               else
                   disp('File generated!')
                   params.y = feature_mat_path.pred_target(:, ti);
                   params.outdir = fullfile(out_dir, process_queue{pr_i}, feature_mat_path.pred_klist{ti}, ['seed_' num2str(seed)]);
               end
               if exist(out_file)
                   a = load(out_file);
                   for scorei = 1:1%size(params.y, 2)
                        a_pfm = load([params.outdir '/PFM_score' num2str(scorei) '_all_folds.mat']);
                        
                        if seed == 1
                                if scorei == 1
                                    eval([process_queue{pr_i} '_pfm_cell = cell(size(params.y,2), 1);']);
                                end
                                eval([process_queue{pr_i} '_pfm_cell{scorei} = a_pfm.PFM_all_folds;']);
                        else
                                eval([process_queue{pr_i} '_pfm_cell{scorei} = [ ' process_queue{pr_i} '_pfm_cell{scorei} a_pfm.PFM_all_folds];']);
                        end
                   end
                   if seed == 1 
                       eval([process_queue{pr_i} '_acc = {a.optimal_acc};']);
                       
                   else
                       eval([process_queue{pr_i} '_acc{end+1} = a.optimal_acc;'])
                   end
               end

       end
   end
for pr_i = 1:length(process_queue)
   eval( [ 'temp_acc = ' process_queue{pr_i} '_acc;']);
   acc_list = [];
   for kki = 1:length(temp_acc)
        acc_list = [acc_list; temp_acc{kki}];
   end
    
   eval(['save(fullfile(out_dir, process_queue{pr_i}, feature_mat_path.pred_klist{ti}, "prediction_ret.mat"),  "' process_queue{pr_i} '_acc", "' process_queue{pr_i} '_pfm_cell", "acc_list");']);
end

end
%    eval(['ret_mat.' target_queue{target_i} '_acc = ' target_queue{target_i} '_acc;']);


% Optimize resolution
% KRR_projects_set = {fullfile(out_dir, '100'), ...
%                     fullfile(out_dir, '200')};

% % Set up parameters for optimizing resolution
% params_opt.nsplit = '';
% params_opt.num_folds = '3';
% params_opt.num_behaviors = '3';
% params_opt.outstem = 'example_targets';
% params_opt.project_set = KRR_projects_set;
% params_opt.out_dir = fullfile(out_dir, 'opt_res');
% params_opt.lambda_set = params.lambda_set;
% params_opt.ker_param = params.ker_param;
% params_opt.sub_fold = params.sub_fold;

% CBIG_GradPar_KRR_optimize_res_wrapper(params_opt);

% rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'predict_phenotypes', 'Kong2023_GradPar'));
% end