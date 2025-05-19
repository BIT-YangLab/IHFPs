# IHFPs

The code of the individualized homologous functional parcellation procedure

## Directory structure  
<pre> IHFPs/  
├── code/ : MATLAB function and script code  
│ ├── IHFPs_config_prepare.m: environment configuration file  
├── example/  
│ ├── IHFPs_procedure.m: example script file for the IHFPs procedure  
| ├── exe_start_gordon_procedure.m: script for individual boundary map and iterative alignment of boundary map  
| ├── CBIG_example_wrapper_gwMRF.m: script for task-constrained group parcellations  
| ├── CBIG_example_wrapper_cmshbm.m: script for cMS-HBM  
| ├── exe_create_functional_IHFPs.m: script for homologous individual parcellations  
| ├── gamlss_fig_script.m: script to use gamlss model to fit developmental trajectories of global mean functional connectivity  
├── resources/  
│ ├── hcd_age_group/ :  age groups of subjects  
├── README.md   </pre>
  
## Dependencies  
gordon gradient map: https://github.com/MidnightScanClub/MSCcodebase,  
gwMRF: https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal,  
cMS-HBM: https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Kong2022_ArealMSHBM,  
gamlss: https://github.com/sunlianglong/BrainChart-FC-Lifespan,  
behavior prediction: https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/predict_phenotypes/Kong2023_GradPar, https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/predict_phenotypes/ChenTam2022_TRBPC,  
MSM: https://github.com/ecr05/MSM_HOCR  
