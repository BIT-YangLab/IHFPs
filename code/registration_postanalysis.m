function registration_postanalysis(out_dir, threshperc, iter_list, work_dir, file_stem)
% average individual gradient maps to population gradient map, create group parcellation as Gordon2016(Option).  
%% input
% out_dir: output data stored in out_dir
% threshperc: threshold of percent using in parcel_creator_cifti
% iter_list: list of iterative surface alignment
% work_dir: 
% file_stem: set output filename
%% outut
% population gradient map


    disp('initial procedure: create input file');
    if ~exist(out_dir)
        mkdir(out_dir);
    end

    % use registration boundary to create parcel
    registration_out_dir = [out_dir '/boundary'];

    model_cifti_file = [ getenv("RESOURCE_DIR") '/model_avg_corrofcorr_allgrad_LR_subcort_smooth2.55_wateredge_avg.dtseries.nii'];
    load([ getenv("RESOURCE_DIR") '/fs_LR_32k_no_medial_wall_index.mat' ]);
    cifti_struct = ft_read_cifti_mod(model_cifti_file);

    boundary_list = zeros(59412,30);

    for iter_i = 1:length(iter_list)
        iter = iter_list(iter_i);
        iter_str=["_iter" num2str(iter)];
        if iter==0
            iter_str="";
        end
        population_mean_boundary_L = cell2mat([registration_out_dir "/population_avg_corrofcorr_allgrad_LR_subcort_L_wateredge_avg" iter_str ".func.gii"]);
        population_mean_boundary_R = cell2mat([registration_out_dir "/population_avg_corrofcorr_allgrad_LR_subcort_R_wateredge_avg" iter_str ".func.gii"]);
        population_mean_boundary_L_gifti = gifti(population_mean_boundary_L);
        population_mean_boundary_R_gifti = gifti(population_mean_boundary_R);
        population_mean_boundary = [population_mean_boundary_L_gifti.cdata; population_mean_boundary_R_gifti.cdata];
        population_mean_boundary = population_mean_boundary(no_medial_wall_index);
        boundary_list(:,iter+1) = population_mean_boundary;
        cifti_struct.data = population_mean_boundary;
        iter_str = ["_iter" num2str(iter)];
        ft_write_cifti_mod(cell2mat([out_dir '/population_avg_corrofcorr_allgrad_LR_subcort_smooth2.55_wateredge_avg_iter' iter_str '.dtseries.nii']), cifti_struct);
        if 0
            parcel_creator_cifti(cell2mat([out_dir '/population_avg_corrofcorr_allgrad_LR_subcort_smooth2.55_wateredge_avg_iter' iter_str '.dtseries.nii']),cell2mat([out_dir file_stem iter_str]), threshperc);
        end
    end

end

