


%% prepare
GORDON_DIR = getenv("GORDON_DIR");
DATASET_DIR = getenv("DATASET_DIR");
HFIP_CODE_DIR = getenv("HFIP_CODE_DIR");
RESOURCE_DIR = getenv("RESOURCE_DIR");
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');


WORK_DIR = [CBIG_CODE_DIR '/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/'];

load([HFIP_CODE_DIR '/resource/fs_LR_32k_no_medial_wall_index.mat']);
load([RESOURCE_DIR '/subject_info.mat']);


%% surface template
targ_mesh = 'fs_LR_32k';

%% input prepare
out_path = ([ getenv('OUT_RET_DIR') '/gwMRF/task_constraind_MRF/' ]);

for age_i = 7:3:19

    output_folder = [ out_path '/HCD_age_related_' num2str(age_i)];
    sub_list_path = [ RESOURCE_DIR  '/hcd_age_group/hcd_sub_qc_new_age' num2str(age_i) '.txt'];
    sub_list = sub_info.src_subject_id(sub_info.group_idx == ((age_i-7)/3+1));
    writecell(sub_list, sub_list_path);
    
    
    sub_num = length(sub_list);
    
    gradient_file = [ RESOURCE_DIR '/model_avg_corrofcorr_allgrad_LR_subcort_smooth2.55_wateredge_avg.dtseries.nii'];
    
    lh_grad_file = [ WORK_DIR '/Code/lib/input/lh_borders_gordon_age_' num2str(age_i) '.mat'];
    rh_grad_file = [ WORK_DIR '/Code/lib/input/rh_borders_gordon_age_' num2str(age_i) '.mat'];
    
    %% create input file
    cmd = ['source ' WORK_DIR ...
        '/examples/example_input/CBIG_gwMRF_create_example_input_fullpaths_1.sh  ' output_folder ...
        '  ' sub_list_path '  0 ' num2str(sub_num-1) ' ' ];
    system(cmd);
    
    
    %% create task-evoked map
    load([ RESOURCE_DIR '/fs_LR_32k_no_medial_wall_index.mat' ]);
    
    task_list = {'EMOTION', 'CARIT', 'GUESSING'};
    cope_list = 1:6;
    task_evoke_map = zeros(59412,18);
    nc = 1;
    for tk_l = 1:length(task_list)
        for cp_i = 1:length(cope_list)
            cope_f = [ DATASET_DIR '/HCD/Package_1188947_hcpdtaskfmriRecom/age_related_' num2str(age_i) '/' task_list{tk_l} '/GrayordinatesStats/cope' num2str(cope_list(cp_i)) '.feat/zstat1.dtseries.nii'];
            task_beta = ft_read_cifti(cope_f);
            dt = task_beta.dtseries(no_medial_wall_index, :);
    
            index =  find(dt >= 5.01);
            task_evoke_map(index, nc) = dt(index);
    
            nc = nc + 1;
        end
    end
    
    save([output_folder '/task_evoke'], 'task_evoke_map');
    
    %% create gradient file
    if contains(targ_mesh, 'fs_LR_32k')
    
        cifti_s = ft_read_cifti_mod(gradient_file);
        
        gradient_label = zeros(64984, 1);
        gradient_label(no_medial_wall_index) = cifti_s.data;
        lh_grad = gradient_label(1:32492);
        rh_grad = gradient_label(32493:end);
        addpath([ WORK_DIR '/Code/lib']);
        [lh_grad_matrix,rh_grad_matrix] = alex_gradient_vertices_to_matrix(lh_grad, rh_grad, targ_mesh);
    
        border_matrix = lh_grad_matrix';
        border = mean(lh_grad_matrix',1);
        save(lh_grad_file,'border_matrix','border');
    
        border_matrix = rh_grad_matrix';
        border = mean(rh_grad_matrix',1);
        save(rh_grad_file,'border_matrix','border');
    end
        
    
    example_input_fullpaths = fullfile(output_folder,'example_input_fullpaths.csv');
    
    code_dir = [WORK_DIR '/Code'];
    cd(code_dir);
    
    %% run procedure
    CBIG_gwMRF_build_data_and_perform_clustering(example_input_fullpaths,output_folder,1,sub_num,200,200,150000,60,1,50000000,15, targ_mesh, lh_grad_file, rh_grad_file);
    
    
end
    



