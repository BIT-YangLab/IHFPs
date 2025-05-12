

%% prepare
compute_orig_homo_rsfc_mean

out_dir = [ getenv('OUT_RET_DIR') '/'];
load([out_dir '/group_network_large_orig_homo.mat']);

ret_info_all = cell(1, 2);

load([ out_dir '/group_17_global_mean_fc.mat');
ret_info_all{1, 1} = fc_parcel_list_all;

load([out_dir '/homologous_global_mean_fc.mat']);
ret_info_all{1, 2} = fc_parcel_list_all;




s_list = {'large', 'homo'};

k_name = 'fc_all';


for si = 1:2
    d_curv = load([out_dir '/dev_curv_init_info.mat']);
    parcel_list = [];
    s_name = sprintf('dev_%s_%s', k_name, s_list{si});
    dst_file_dir = [out_dir '/prediction_ret/' s_name];
    if ~exist(dst_file_dir)
        mkdir(dst_file_dir);
    end
    dst_mat_f = [out_dir '/curv_info/' s_name '.mat'];

    ret_cell_list = ret_info_all{si}';
    nc = 1;
    for pi = 1:size(ret_cell_list, 1)
        tmp_s = d_curv.rsfc_mean;
        tmp_s(:, 1) = ret_cell_list(pi, :);
        idx = union(find(tmp_s(:, 1) == 0), find(isnan(tmp_s(:, 1))));

        tmp_s(idx, :) = [];
        if isempty(tmp_s)
            continue;
        end
        parcel_list(end+1) = pi;
        eval(['d_curv.rs_info_cnt' num2str(nc) ' = tmp_s;']);
        nc = nc + 1;
        
    end
    d_curv.parcel_list = parcel_list;
    save(dst_mat_f, 'd_curv');



end


d_curv = load('/home/jinlong/VDisk2/Jinlong/2024_homologous_parcellation_result/HFIP_result/HFIP_ret/dev_curv_init_info.mat');
parcel_list = [];
s_name = sprintf('dev_%s_%s_mean_network_weight', k_name, s_list{si});
dst_file_dir = ['/home/jinlong/VDisk2/Jinlong/2024_homologous_parcellation_result/HFIP_result/prediction_ret/' s_name];
if ~exist(dst_file_dir)
    mkdir(dst_file_dir);
end
dst_mat_f = ['/home/jinlong/VDisk2/Jinlong/2024_homologous_parcellation_result/HFIP_result/HFIP_ret/curv_info/' s_name '.mat'];

all_sub_info = zeros(591, 17);
%         eval(['ret_cell_tmp = ' s_list{si} '_ret_cell;']);
ret_cell_list = ret_info_all{1, si};




for subi = 1:591
    ta = ret_cell_list(subi, :);
    ta(isnan(ta)) = 0;
    tb = all_sub_assign_prob{subi};
    tc = ta * tb;
    tm = sum(tb, 1);
    tm(tm == 0) = 1;
    tc = tc ./ tm;
    all_sub_info(subi, :) = tc;
end


for neti = 1:17
    tmp_s = d_curv.rsfc_mean;
    tmp_s(:, 1) = all_sub_info(:, neti);
    idx = union(find(isnan(tmp_s(:, 1))), find(tmp_s(:, 1) == 0));
    tmp_s(idx, :) = [];

    if isempty(tmp_s)
        continue;
    end

    eval(['d_curv.rs_info_cnt' num2str(neti) ' = tmp_s;']);
end
save(dst_mat_f, 'd_curv');



%% Rscript GAMLSS model
% fit areal-level and network-level normative trajectories
system(['Rscript ' getenv('HFIP_CODE_DIR') '/create_gamlss_info.R']);

% variability of areal-level and network-level normative trajectoris
compute_trajectories_variability('dev_fc_all.mat', 'dev_fc_all', 'dev_fc_all_mean_network_weight');

% growth rate of network-level normative trajectoris
system(['../code/bootstrapping_script.sh']);



