clc;clear;

% example script of Homologous individual parcellations

% step1: configure
addpath('/home/jinlong/VDisk1/Jinlong/2024_IHFPS/code');
IHFPs_config_prepare

% step2: individual boundary map and iterative alignment of boundary map
exe_start_gordon_procedure

% step2: task-constrained group parcellations
CBIG_example_wrapper_gwMRF

% step4: fine-grained individual parcellations
CBIG_example_wrapper_cmshbm

% step5: homologous individual parcellations
exe_create_functional_IHFPs

% step6: mapping parcellations to large-scale networks
out_dir = [ getenv('OUT_RET_DIR') '/network_level_ret'];
load([ getenv('OUT_RET_DIR') '/IHFP_ret/IHFP_sess1.mat'] );

group_fine2network_profile(out_dir, group_fine_label);

group2ind_TM_network_profile(out_dir, group_fine_label, new_label_list_all, all_sub_list);

homologous_ind_match(out_dir, group_fine_label, new_label_list_all, all_sub_list);


% step7: using gamlss model to fit developmental trajectories of global mean fc
gamlss_fig_script



% step8: predict behavior measures using RSFC within and between ROI  
out_dir = '/nd_disk3/guoyuan/Jinlong/2024IHFPs/prediction_ret/pred_rsfc_between' ; 
prediction_krr_behavior(out_dir, {'pred_rsfc_homo', 'pred_rsfc_orig'}, [out_dir '/pred_rsfc_between.mat'])
pfm_significant_matrix_create
pfm_significant_matrix_create_permutation
pfm_rsfc_map_create





