

function create_population_boundary(out_dir, sub_path, iter_list, boundary_path)
% create population boundary map
%% input:
% out_dir:
% sub_path: path of subjects' list
% iter_list: list of iterative surface alignment
% boundary_path: source path of individual boundary maps
%% output:
% population boundary maps in 'out_dir/boundary/'

    out_file_dir=[out_dir '/boundary'];
    if ~exist(out_file_dir)
        mkdir(out_file_dir);
    end
    gradsname = 'avg_corrofcorr_allgrad_LR_subcort';


    L_population_mean_func=zeros(32492,60);
    R_population_mean_func=zeros(32492,60);

    sublist=textread(sub_path,'%s');

    for iter_idx = 1:length(iter_list)
        iter_str = ['_iter' num2str(iter_list(iter_idx))];
        if iter_list(iter_idx) == 0
            iter_str = '';
        end
        disp(['process iter:' num2str(iter_idx)]);
        for i=1:length(sublist)
            disp(['process sub#' num2str(i)]);
            boundary_file_L_path=cell2mat([boundary_path '/boundary/' sublist(i) '_avg_corrofcorr_allgrad_LR_subcort_L_wateredge_avg_mod' iter_str '.func.gii']);
            boundary_file_R_path=cell2mat([boundary_path '/boundary/' sublist(i) '_avg_corrofcorr_allgrad_LR_subcort_R_wateredge_avg_mod' iter_str '.func.gii']);
            boundary_file_L=gifti(boundary_file_L_path);
            boundary_file_R=gifti(boundary_file_R_path);
            boundary_L=boundary_file_L.cdata;
            boundary_R=boundary_file_R.cdata;
            L_population_mean_func(:,i)=boundary_L;
            R_population_mean_func(:,i)=boundary_R;
        end
        
        L_population_mean_func=mean(L_population_mean_func,2);
        R_population_mean_func=mean(R_population_mean_func,2);
        L_func_file=([out_file_dir '/population_' gradsname '_L_wateredge_avg' iter_str '.func.gii']);
        R_func_file=([out_file_dir '/population_' gradsname '_R_wateredge_avg' iter_str '.func.gii']);
        save(gifti(single(L_population_mean_func)), L_func_file);
        save(gifti(single(R_population_mean_func)), R_func_file);
        system(['wb_command -set-structure ' L_func_file ' CORTEX_LEFT']);
        system(['wb_command -set-structure ' R_func_file ' CORTEX_RIGHT']);

    end




end








