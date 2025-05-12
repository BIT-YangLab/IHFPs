%% set path
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
OUT_RET_DIR = [ getenv('OUT_RET_DIR') '/Kong2022_ind/' ];
work_dir = [CBIG_CODE_DIR '/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/'];
seed_mesh='fs_LR_900';
targ_mesh='fs_LR_32k';
split_flag = '0';
sess_list = 1:1;
num_sess = 1;
num_clusters = 400;
beta = '100';
num_estimate_s = 10;
w='100';
c='50';

load([getenv("RESOURCE_DIR") '/subject_info.mat']);

%% process
for age_i = 7:3:19

    out_dir = [ OUT_RET_DIR '/Kong2022_ret/HCD_age_group_' num2str(age_i)];

    cur_dir = pwd;
    if ~exist(out_dir)
        mkdir(out_dir);
    else
        rmdir(out_dir, 's');
    end

    sub_path = [ getenv("RESOURCE_DIR") '/hcd_age_group/hcd_sub_qc_new_age' num2str(age_i) '.txt'];
    sub_list = sub_info.src_subject_id(sub_info.group_idx == ((age_i-7)/3+1));
    writecell(sub_list, sub_path);

     
    sub_num = length(sub_list);
    

    

    cmd = [work_dir '/examples/CBIG_create_input_data_mod.sh   ' out_dir   '   '  sub_path ];
    system(cmd);

    %% generate profile
    cd ([work_dir '/step1_generate_profiles_and_ini_params']);
    step_out_dir = [out_dir '/generate_profiles_and_ini_params' ];
    for i = 1:sub_num
        sub = sub_list{i,1};
        for sess = 1:num_sess
            CBIG_ArealMSHBM_generate_profiles(seed_mesh,targ_mesh,step_out_dir,sub, num2str(sess),split_flag);
        end
    end

    CBIG_ArealMSHBM_avg_profiles(seed_mesh, targ_mesh, step_out_dir, ceil(sub_num/3), num_sess, sub_list);
   
    a = load([ getenv('OUT_RET_DIR') '/gwMRF/task_constraind_MRF/HCD_age_related_' num2str(age_i) '/clustering/Graph_Cut_faster__grad_prior_gordon_cluster_20_datacost_1_smoothcost_1500_iterations_60_seed_1.mat']);
    lh_label = a.results.lh_label';
    rh_label = a.results.rh_label';

    CBIG_ArealMSHBM_generate_ini_params(seed_mesh, targ_mesh, lh_label, rh_label, step_out_dir);


    cmd = [ work_dir '/examples/CBIG_create_input_data_step1.sh   ' out_dir   '   '  sub_path ' ' num2str(num_estimate_s)];
    system(cmd);

    CBIG_ArealMSHBM_generate_radius_mask(lh_label, rh_label, targ_mesh, '30', step_out_dir);

    %% estimate prior
    step_out_dir = [out_dir '/estimate_group_priors' ];
    cd ([work_dir '/step2_estimate_priors']);

    cmd = [ work_dir '/examples/CBIG_create_input_data_step2.sh   ' out_dir   '   '  sub_path];
    system(cmd);

    tmp_dir = [out_dir '/estimate_group_priors/tmp_results1'];
    
    fid = fopen([ out_dir '/estimate_script1.sh'], 'w');
    for i=1:num_estimate_s
        disp(i);
        fprintf(fid, '%s\n', ['matlab -nojvm -nodesktop -nosplash -r "pause(10); addpath('    char("[ work_dir '/step2_estimate_priors']);  addpath([ work_dir '/examples']); CBIG_config_path;  CBIG_ArealMSHBM_cMSHBM_estimate_group_priors_child(")  char("'") num2str(i) char("','") targ_mesh char("','") tmp_dir char("'") '); exit;"   &']);
    end
    fprintf(fid, '%s\n' , ['matlab -nojvm -nodesktop -nosplash -r " addpath('    char("[ work_dir '/step2_estimate_priors']); addpath([ work_dir '/examples']); CBIG_config_path; CBIG_ArealMSHBM_cMSHBM_estimate_group_priors_parent(")  char("'") step_out_dir char("','") targ_mesh char("','") num2str(num_estimate_s) char("','") num2str(1) char("','") '400' char("','") beta char("','") tmp_dir char("','") '50' char("'") '); exit;"  ']);
    fclose(fid);
    system(['chmod a+x  ' out_dir '/estimate_script1.sh']);
    system([ out_dir '/estimate_script1.sh']);

    %% generate individual parcellations
    step_out_dir = [out_dir '/generate_individual_parcellations' ];
    cd ([work_dir '/step3_generate_ind_parcellations']);
    cmd = [ work_dir '/examples/CBIG_create_input_data_step3.sh   ' out_dir   '   '  sub_path  '  cMSHBM ' num2str(beta)];
    system(cmd);


    

    for st = 1:num_sess
        for sub_i = 1:sub_num
            [lh_labels, rh_labels] = CBIG_ArealMSHBM_cMSHBM_generate_individual_parcellation( ...
            step_out_dir, targ_mesh, num2str(1), num2str(num_clusters), num2str(sub_i), (w), (c), beta, 'validation_set', st);
            if sub_i == 1
                label_list = [lh_labels;rh_labels];
            else
                label_tmp = [lh_labels;rh_labels];
                label_list = [label_list label_tmp];
            end
        end
        save([step_out_dir '/label_list_' num2str(st) '.mat'], 'label_list');
    end
    a = 1;
end