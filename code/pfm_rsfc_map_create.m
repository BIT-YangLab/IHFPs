
format_rsfc_behv = [ getenv('OUT_RET_DIR') '/prediction_ret/pred_rsfc_between/%s/%s/prediction_ret.mat'];

load([ getenv('OUT_RET_DIR') '/prediction_initial_info_new.mat']);
load([ getenv('OUT_RET_DIR') '/group_network_large_orig_homo.mat']);

sk_list = {'pred_rsfc_orig', 'pred_rsfc_homo'};
parcel_n_list = [400, 400];
pro_name_list = {'orig', 'homo'};

% Yeo 17 network
label_name = {'Auditory','Dorsal Attention A','Control A', 'Somatomotor A', 'Salience/VenAttn B', ...
    'Default B', 'Default C', 'Visual C', 'visual A', 'Dorsal Attention B', 'Temoral Pariental', 'Control B','Visual B',...
    'Control C','Default A','Salience/VenAttn A', 'Somatomotor B'};
pfm_resort_list = [11 15 6 7 3 12 14 1 16 5 2 10 4 17 9 13 8];
pfm_resort_list_l = [1 4 7 8 10 12 14 16 ];

perm_num = 2000;
perm_cnt = 40;


match_plist_orig = [];
match_plist_homo = [];


net_match_plist_orig = [];
net_match_plist_homo = [];


split_list_orig = [];
split_list_homo = [];




nc = 0;
for neti = 1:17
    match_plist_orig = [match_plist_orig; orig_system_list{pfm_resort_list(neti)}];
    net_match_plist_orig = [net_match_plist_orig; ones(length(orig_system_list{pfm_resort_list(neti)}), 1).*pfm_resort_list(neti)];
    split_list_orig(end+1) = nc + length(orig_system_list{pfm_resort_list(neti)});
    nc = nc + length(orig_system_list{pfm_resort_list(neti)});
end

nc = 0;
for neti = 1:17
    match_plist_homo= [match_plist_homo; homo_system_list{pfm_resort_list(neti)}'];
    net_match_plist_homo= [net_match_plist_homo; ones(length(homo_system_list{pfm_resort_list(neti)}), 1).*pfm_resort_list(neti)];
    split_list_homo(end+1) = nc + length(homo_system_list{pfm_resort_list(neti)});
    nc = nc + length(homo_system_list{pfm_resort_list(neti)});
end

match_pcell = {match_plist_orig, match_plist_homo};
net_match_pcell = {net_match_plist_orig, net_match_plist_homo};
split_pcell = {split_list_orig, split_list_homo};
system_list_cell = {orig_system_list, homo_system_list};
pfm_ret = cell(8, 2);

load( [ getenv('OUT_RET_DIR') '/HFIP_result/prediction_ret/pred_rsfc_between.mat']);

for ki = 1:8

    for si = 1:2
        rsfc_b_f  = sprintf(format_rsfc_behv, sk_list{si}, pred_klist{ki});
        rsfc_b_s = load(rsfc_b_f);
        eval(['pfm = nanmean(rsfc_b_s.' sk_list{si} '_pfm_cell{1}, 2);']);

        eval(['resample_idx = resample_idx_' pro_name_list{pi} ';']);

        net_corr_sign_f = zeros(17, 17);
        pfm1 = zeros(400, 400);
        pfm1(resample_idx) = pfm;
        
        std1 = std(pfm);
        pfm1 = pfm ./ std1;
        pfm1 = reshape(pfm1, parcel_n_list(si), parcel_n_list(si));


        pfm2 = pfm1(match_pcell{si}, match_pcell{si});

        pfm_ret{ki, si} = pfm2;
    end


end



colors = [
    0 1 1;   % Green
    0 0 1;   % Blue
    0 0 0;   % Black
    1 0 0;   % Red
    1 1 0    % Yellow
    ];


nColors = 256;
nSegments = size(colors, 1) - 1;


customColormap = [];


for i = 1:nSegments
    startColor = colors(i, :);
    endColor = colors(i + 1, :);
    segmentColors = interp1([0 1], [startColor; endColor], linspace(0, 1, nColors/nSegments));
    customColormap = [customColormap; segmentColors]; 
end

for ki = 1:8
    for si = 1:3
        fig = figure;
        set(fig, 'Position', [100, 100, 800, 600]);
        colormap(customColormap);
        imagesc(pfm_ret{ki, si});
        colorbar;
        xticks([]);
        yticks([]);
        max_abs_value = max(abs(pfm_ret{ki, si}(:)));max_abs_value = 2.5;
        caxis([-max_abs_value, max_abs_value]);
        for gi = 1:length(split_pcell{si})
            hold on;
            if nnz(pfm_resort_list_l == gi)
                line([0.5 max(split_pcell{si}) + 0.5], [split_pcell{si}(gi)+0.5 split_pcell{si}(gi)+0.5], 'Color', 'white', 'LineWidth', 0.5);
                line([split_pcell{si}(gi)+0.5 split_pcell{si}(gi)+0.5], [0.5 max(split_pcell{si})+0.5], 'Color', 'white', 'LineWidth', 0.5);

            else
                line([0.5 max(split_pcell{si})+0.5], [split_pcell{si}(gi)+0.5 split_pcell{si}(gi)+0.5], 'Color', 'white', 'LineWidth', 0.001);
                line([split_pcell{si}(gi)+0.5 split_pcell{si}(gi)+0.5], [0.5 max(split_pcell{si})+0.5], 'Color', 'white', 'LineWidth', 0.001);
            end


       
        end


        hold off
        print([getenv('OUT_RET_DIR') '/pic_dir/pfm_rsfc/pfm_' sk_list{si} '_' pred_klist{ki} '.eps'],'-depsc','-r600');% save eps with 600dpi

        pause(5)
    end

end
