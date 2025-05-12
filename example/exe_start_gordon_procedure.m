


%% prepare
GORDON_DIR = getenv("GORDON_CODE_DIR");
DATASET_DIR = getenv("DATASET_DIR");
HFIP_CODE_DIR = getenv("HFIP_CODE_DIR");
RESOURCE_DIR = getenv("RESOURCE_DIR");

addpath([ GORDON_DIR '/cifti-matlab']);
addpath([ GORDON_DIR '/cifti-matlab-WashU-gradient']);


%% create age-independent boundary map
% example test
max_iter = 15;
out_dir = [  getenv("OUT_RET_DIR") '/Gordon2016_mod/HCD_591qc_age_independent/']; 
cohort_file = [out_dir '/cohort.txt'];
tmask_file = [out_dir '/tmasklist.txt'];

load([RESOURCE_DIR '/subject_info.mat']);
subject_list = sub_info.src_subject_id;
sub_path = [ RESOURCE_DIR '/hcd_age_group/hcd_sub_qc_new_age_independent.txt'];
writecell(subject_list, sub_path);

%% data initial
disp('create initial file...');
initial_procedure(sub_path, out_dir, cohort_file, tmask_file);

%% gradient map
disp('create individual gradient file');
surface_parcellation(cohort_file,tmask_file,50,0,out_dir, work_dir);
create_population_boundary(out_dir, sub_path, 0, out_dir);
system(['cd  ' HFIP_CODE_DIR '/surface_registration; ./start_registration.sh ' HFIP_CODE_DIR '/surface_registration  ' sub_path ' '  out_dir  ' ' num2str(max_iter) ' ' out_dir ' ; cd ' work_dir]);



for age_i = 1:5
    % subject list file
    disp('create sub list');
    
    sub_path = [ RESOURCE_DIR '/hcd_age_group/hcd_sub_qc_new_age' num2str(age_i ) '.txt'];
    subject_list = sub_info.src_subject_id(sub_info.group_idx == age_i);
    writecell(subject_list, sub_path);

    work_dir = [ GORDON_DIR '/'];

    % for infomap running ( infomap can't run automatically for matlab's gcc env version, so using shell-gcc)
    fid = fopen([work_dir '/Parcels/run_script.sh'], 'w');

    % out directory
    out_dir = [  getenv("OUT_RET_DIR") '/Gordon2016_mod/HCD_591qc_age_' num2str(age_i) '/'];
    boundary_dir = [getenv("OUT_RET_DIR") '/Gordon2016_mod/HCD_591qc_age_independent/'];

    % gradient file path
    if ~exist([out_dir '/boundary'])
        mkdir ([out_dir '/boundary']);
    end


    file_stem = '/Parcels_population_';
    threshperc = 0.60;

    % example test
    max_iter = 15; 
    cohort_file = [out_dir '/cohort.txt'];
    tmask_file = [out_dir '/tmasklist.txt'];

    %% data initial
    disp('create initial file...');
    initial_procedure(sub_path, out_dir, cohort_file, tmask_file);

    %% gradient map
    disp('create individual gradient file');
    surface_parcellation(cohort_file,tmask_file,50,0,out_dir, work_dir);

    %% population boundary
    disp('create population boundary...');
    create_population_boundary(out_dir, sub_path, max_iter, boundary_dir);

    %% parcel create 
    disp('create parcellation...');
    registration_postanalysis(out_dir, threshperc, max_iter, work_dir, file_stem);
    
    

end