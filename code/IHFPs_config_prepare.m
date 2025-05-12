%% prepare running configure
% 1. HFP path, CBIG path, Gordon path, data path, et.
% 2. software path: connectome workbench, fsl, mirtk, msm, et
% 3. 
% 4. 


setenv('HFIP_CODE_DIR', '/home/jinlong/VDisk1/Jinlong/2024_IHFPs');
setenv('CBIG_CODE_DIR', '/home/jinlong/VDisk1/Jinlong/external_package/parcellation_method/CBIG');
setenv('GORDON_CODE_DIR', '/home/jinlong/VDisk1/Jinlong/external_package/parcellation_method/Gordon2016_mod');
setenv('DATASET_DIR', '/home/jinlong/VDisk2/Jinlong/data');
setenv('WB_DIR', '/opt/workbench/bin_linux64/wb_command');
setenv('FSL_DIR', '/usr/local/fsl/share/fsl/bin/fsl');
setenv('FREESURFER_HOME', '/usr/local/freesurfer/6.0.1');
setenv('OUT_RET_DIR', '/home/jinlong/VDisk2/Jinlong/2024_homologous_parcellation_result/HFIP_result');
setenv('RESOURCE_DIR', '/home/jinlong/VDisk1/Jinlong/2024_IHFPs/1resources');

setenv("GORDON_DIR", '/home/jinlong/VDisk1/Jinlong/external_package/parcellation_method/Gordon2016_mod');
setenv("DATASET_DIR", '/home/jinlong/VDisk2/Jinlong/data/');


HFIP_CODE_DIR = getenv('HFIP_CODE_DIR');
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
GORDON_CODE_DIR = getenv('GORDON_CODE_DIR');

addpath([ HFIP_CODE_DIR, '/code']);
addpath([ HFIP_CODE_DIR, '/conn_freq_hcp']);
addpath([ HFIP_CODE_DIR, '/postproc']);
addpath([ HFIP_CODE_DIR, '/external_mod']);


addpath([ CBIG_CODE_DIR, '/utilities/matlab/speedup_gradients']);
addpath([ CBIG_CODE_DIR, '/utilities/matlab/fslr_matlab']);
addpath([ CBIG_CODE_DIR, '/external_packages']);
addpath([ CBIG_CODE_DIR, '/external_packages/SD/SDv1.5.1-svn593/BasicTools/']);
addpath([ CBIG_CODE_DIR, '/utilities/matlab/speedup_gradients/utilities/']);
addpath([ CBIG_CODE_DIR, '/external_packages/SD/SDv1.5.1-svn593/kd_tree/']);
addpath([ CBIG_CODE_DIR, '/external_packages/matlab/default_packages/others/']);
addpath([ CBIG_CODE_DIR, '/utilities/matlab/FC/']);
addpath([ CBIG_CODE_DIR, '/utilities/matlab/stats/']);
addpath([ CBIG_CODE_DIR, '/utilities/matlab/surf/']);
addpath([ CBIG_CODE_DIR, '/external_packages/matlab/default_packages/graph_cut/matlab/']);
addpath([ CBIG_CODE_DIR, '/external_packages/matlab/default_packages/matlab_bgl/']);
addpath([ CBIG_CODE_DIR, '/external_packages/matlab/default_packages/mtimesx_20110223/']);
addpath([ CBIG_CODE_DIR, '/utilities/matlab/parcellation/']);
addpath([ CBIG_CODE_DIR, '/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Code/lib']);
addpath([ CBIG_CODE_DIR, '/external_packages/matlab/default_packages/DSP/']);
addpath([ CBIG_CODE_DIR, '/utilities/matlab/DSP/']);
addpath([ CBIG_CODE_DIR, '/utilities']);
addpath([ CBIG_CODE_DIR, '/utilities/matlab/predictive_models/PFM']);
addpath([ CBIG_CODE_DIR, '/stable_projects/predict_phenotypes/Kong2023_GradPar']);
addpath([ CBIG_CODE_DIR, '/utilities/matlab/predictive_models/utilities']);
addpath([ CBIG_CODE_DIR, '/utilities/matlab/utilities']);
addpath([ CBIG_CODE_DIR, '/utilities/matlab/predictive_models/KernelRidgeRegression']);
addpath([ CBIG_CODE_DIR, '/utilities/matlab/predictive_models/LinearRidgeRegression']);
addpath([ CBIG_CODE_DIR, '/utilities/matlab/stats']);
addpath([getenv('CBIG_CODE_DIR'),  '/external_packages/matlab/default_packages/WashU_gradients/']);

addpath([ GORDON_CODE_DIR, '/code/']);
addpath([ GORDON_CODE_DIR, '/Infomap_wrapper/']);

addpath([ HFIP_CODE_DIR '/1code/compare_matrices_to_assign_networks-main/']);
addpath([ HFIP_CODE_DIR '/1code/plotting-tools-master/custom_hist/']);
addpath([ HFIP_CODE_DIR '/1code/GRETNA/NetFunctions/'])

