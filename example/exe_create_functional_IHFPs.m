

%% configuration
dataset = 'HCD';
mesh = 'fs_LR_32k';
num_session = 1;


sub_list_all_file = [ getenv("RESOURCE_DIR") '/hcd_age_group/hcd_sub_qc_new.txt'];
sub_list_age_file = [ getenv('RESOURCE_DIR') '/hcd_age_group/hcd_sub_qc_new_age%d.txt'];


process_queue = { '7_9' '10_12' '13_15' '16_18' '19_21'};
process_queue_id = [7, 10, 13, 16, 19];



group_label_path = [ getenv('OUT_RET_DIR') '/gwMRF/task_constraind_MRF/clustering/Graph_Cut_faster__grad_prior_gordon_cluster_20_datacost_1_smoothcost_1500_iterations_60_seed_1.mat'];

out_dir = [ getenv('OUT_RET_DIR') '/HIP'];
if ~exist(out_dir)
    mkdir(out_dir)
end

no_medial_wall_index = find_cortical_label_index(mesh);

flg_optimize_ind_atlas = 0;
flg_create_group_ref = 0;
flg_area_match_homo_atlas = 1;
flg_analysis = 0;

load([getenv("RESOURCE_DIR") '/subject_info.mat']);


%% Kong2022 individual parcellation --> fill empty area and remove noise vertices
disp('fill empty area and remove noise vertices');
if flg_optimize_ind_atlas == 1
    new_label_list_all = cell(length(process_queue), num_session);
    for sess = 1:num_session
        eval(['new_label_list_all_sess' num2str(sess) ' = {};' ]);
        disp([' sess ' num2str(sess) ]);
        for p_i = 1:length(process_queue)
            disp(['process ' process_queue{p_i}]);
            load([ getenv('OUT_RET_DIR')  '/Kong2022_ind/HCD_age_related_' num2str(process_queue_id(p_i)) '/generate_individual_parcellations/label_list_' num2str(p_i) '.mat']); % label_list
            sub_list_file = sprintf(sub_list_age_file, process_queue_id(p_i));
            grp_idx = find(sub_info.group_idx == p_i);
            sub_list = sub_info.src_subject_id(grp_idx);
            writecell(sub_list, sub_list_file);
            new_label_list = label_list;
            
            for sub_i = 1:size(label_list, 2)
                disp([num2str(sub_i) '/' num2str(size(label_list,2)) '  ' sub_list{sub_i, 1}]);
                if num_session == 1
                    cifti_file = sub_info.file_all{grp_idx(sub_i)};
                else
                    cifti_file = sub_info.file_sess2{grp_idx(sub_i)}{sess};
                end
                [new_label_list(:, sub_i), k_flag] = erase_noise_vertice(new_label_list(:, sub_i), mesh);
                new_label_list(:, sub_i) = fill_empty_area(label_list(:, sub_i), {cifti_file}, mesh);

            end
            new_label_list_all{p_i, sess} = new_label_list;
        end
    end
    save([ getenv("OUT_RET_DIR") '/IHFP_ret/new_label_list_all_sessall.mat'], 'group_label_path', 'new_label_list_all');
else
    load([ getenv("OUT_RET_DIR") '/IHFP_ret/new_label_list_all_sessall.mat']);
end

% -> new_label_list_all


%% create group-level ref label

group_result = load(group_label_path);
lh_label = group_result.results.lh_label';
rh_label = group_result.results.rh_label';
rh_label(rh_label > 0) = rh_label(rh_label > 0) + max(lh_label);
age_independent_group_label = [lh_label; rh_label];



% remove small parcel (vertices less than 4) 
parcel_list = unique(age_independent_group_label);
for i = 1:length(parcel_list)
    index = find(age_independent_group_label == parcel_list(i));
    if length(index) < 4
        age_independent_group_label(index) = 0;
    end
end

sub_list = sub_info.src_subject_id;
writecell(sub_list, sub_list_all_file);
cifti_all_path = {};

if flg_create_group_ref == 1

    if nnz(age_independent_group_label(no_medial_wall_index) == 0)
        for sess = 1:num_session
            for sub_i = 1:length(sub_list)
                if num_session == 1
                    cifti_file = sub_info.file_all{(sub_i)};
                else
                    cifti_file = sub_info.file_sess2{(sub_i)}{sess};
                end
                
                cifti_all_path{end+1} = cifti_file;
            end
        end
        group_label = fill_empty_area(age_independent_group_label, cifti_all_path, mesh);
    else
        group_label = age_independent_group_label;
    end
    group_label = erase_noise_vertice(group_label, mesh);

    group_ref_label = create_ref_groupROI('./', group_label, new_label_list_all);
else
    load([getenv('OUT_RET_DIR') '/IHFP_ret/IHFP_group_ref_label.mat'])

end
% -> group_ref_label

%% match homologous network
disp('match homologous neiwork');
if flg_area_match_homo_atlas == 1
    for sess = 1:num_session
        for p_i = 1:length(process_queue)
            disp(['process ' process_queue{p_i} ]);
            sub_list_file = sprintf(sub_list_age_file, process_queue_id(p_i));
            grp_idx = find(sub_info.group_idx == p_i);
            sub_list = sub_info.src_subject_id(grp_idx);
            writecell(sub_list, sub_list_file);
            label_list = new_label_list_all{sess, p_i};
            new_label_list = label_list;

            for sub_i = 1:size(label_list, 2)
                disp(['sub ' num2str(sub_i) '/' num2str(length(sub_list)) ' ' sub_list{sub_i}]);
                if num_session == 1
                    cifti_file = sub_info.file_all{(sub_i)};
                else
                    cifti_file = sub_info.file_sess2{(sub_i)}{sess};
                end
                
                ind_new_label =  match_homologous_network(group_ref_label, label_list(:, sub_i), {cifti_file}, 1);
                
                ind_new_label = fill_empty_area(ind_new_label, {cifti_file}, mesh);
                ind_new_label = erase_noise_vertice(ind_new_label, mesh);
                new_label_list(:, sub_i) = ind_new_label;
            end
            new_label_list_all{p_i} = new_label_list;
            if p_i == 1
                all_label_list = new_label_list;
                all_sub_list = sub_list;
            else
                all_label_list = [all_label_list new_label_list];
                all_sub_list = [all_sub_list; sub_list];
            end
        end

        %% find stable network across all subjects
        all_label_list_no_filter = all_label_list;
        stable_net = find_stable_homogeneity_patch(all_label_list);
        for sub_i = 1:size(all_label_list, 2)
            ind_label = zeros(size(all_label_list(:, sub_i)));
            for j = 1:length(stable_net)
                index = find(all_label_list(:, sub_i) == stable_net(j));
                ind_label(index) = all_label_list(index, sub_i);
            end
            all_label_list(:, sub_i) = ind_label;
        end
        sess_struct.all_label_list = all_label_list;
        sess_struct.new_label_list_all = new_label_list_all;
        sess_struct.all_label_list_no_filter = all_label_list_no_filter;
        eval(['sess' num2str(sess) '_struct = sess_struct;' ]);

        save([ getenv("OUT_RET_DIR")  '/IHFP_ret/IHFP_sess' num2str(sess) '.mat'], 'age_independent_group_label', 'all_label_list', 'all_sub_list', ...
            'group_label', 'group_ref_label', 'new_label_list_all', 'no_medial_wall_index', 'stable_net', 'all_label_list_no_filter', 'sess_struct');
    end

else
    for sess = 1:num_session
        eval(['sess' num2str(sess) '_struct = load([ getenv("OUT_RET_DIR")  ' char("'") '/IHFP_ret/IHFP_sess' num2str(sess) '.mat' char("'") ']);']);
    end
end








