clc;clear;

src_name = 'dev_fc_homogeneity_homo';
load('/home/jinlong/VDisk2/Jinlong/resources/net17_colorbar.mat');
load(['/home/jinlong/VDisk2/Jinlong/2024_homologous_parcellation_result/HFIP_result/HFIP_ret/curv_info/' src_name '.mat']);
load('/home/jinlong/VDisk2/Jinlong/2024_homologous_parcellation_result/HFIP_result/HFIP_ret/group_network_large_orig_homo.mat');
load('/home/jinlong/VDisk2/Jinlong/2024_homologous_parcellation_result/HFIP_result/HFIP_ret/homologous_large_network_parcel_match_17_new2.mat');
format_dst_dir = '/home/jinlong/VDisk2/Jinlong/2024_homologous_parcellation_result/HFIP_result/prediction_ret/%s/gamlss_fit_allData_rs_global_fc_%d.mat';

label_name = {'Auditory','Dorsal Attention A','Control A', 'Somatomotor A', 'Salience/VenAttn B', ...
                     'Default B', 'Default C', 'Visual C', 'visual A', 'Dorsal Attention B', 'Temoral Pariental', 'Control B','Visual B',...
                     'Control C','Default A','Salience/VenAttn A', 'Somatomotor B'};
pfm_resort_list = [11 15 6 7 3 12 14 1 16 5 2 10 4 17 9 13 8];
pfm_resort_list_l = [1 4 7 8 10 12 14 16 ];

parcel_list = d_curv.parcel_list;

net17_colorbar = net17_colorbar(pfm_resort_list, :);

parcel_n_list = [34, 400, 400];



match_plist_large = zeros(34, 1);
match_plist_orig = [];
match_plist_homo = [];

net_match_plist_large = zeros(34, 1);
net_match_plist_orig = [];
net_match_plist_homo = [];

split_list_large = [];
split_list_orig = [];
split_list_homo = [];




nc = 0;
for neti = 1:17
    match_plist_homo= [match_plist_homo; homo_system_list{pfm_resort_list(neti)}'];
    net_match_plist_homo= [net_match_plist_homo; ones(length(homo_system_list{pfm_resort_list(neti)}), 1).*pfm_resort_list(neti)];
    split_list_homo(end+1) = nc + length(homo_system_list{pfm_resort_list(neti)});
    nc = nc + length(homo_system_list{pfm_resort_list(neti)});
end

match_pcell = {match_plist_large, match_plist_orig, match_plist_homo};
split_pcell = {split_list_large, split_list_orig, split_list_homo};
pfm_ret = cell(8, 3);
pfm_sa_ret = cell(2, 3);

        ad1 = homo_match_ret(:, 2);
[mv, mi] = sort(ad1);
new_match = homo_match_ret(mi, :);

corrmap = paircorr_mod(all_cen_y_list);

corr1(d_curv.parcel_list, d_curv.parcel_list) = corrmap;

% tmp = sort(reshape(corr1, 1, []), 'descend');
% ths_r = tmp(round(length(tmp)*0.1));
% ths_l = tmp(round(length(tmp)*0.9));
% flg_r = corr1 < ths_r;
% flg_l = corr1 > ths_l;
% corr1(boolean(flg_r&flg_l)) = 0;

        corr2 = corr1(match_plist_homo, match_plist_homo);

 



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
    customColormap = [customColormap; segmentColors]; % ½«Éú³ÉµÄÑÕÉ«Ìí¼Óµ½×Ô¶¨ÒåÑÕÉ«Ó³ÉäÖÐ
end
si = 3;
        fig = figure;
        set(fig, 'Position', [100, 100, 800, 600]);
        colormap(customColormap);
        imagesc(corr2);
        colorbar;
        xticks([]); 
        yticks([]); 
        for gi = 1:length(split_pcell{si})
            hold on;
            if nnz(pfm_resort_list_l == gi)
                line([0.5 max(split_pcell{si}) + 0.5], [split_pcell{si}(gi)+0.5 split_pcell{si}(gi)+0.5], 'Color', 'white', 'LineWidth', 0.5);
                line([split_pcell{si}(gi)+0.5 split_pcell{si}(gi)+0.5], [0.5 max(split_pcell{si})+0.5], 'Color', 'white', 'LineWidth', 0.5);

            else
                line([0.5 max(split_pcell{si})+0.5], [split_pcell{si}(gi)+0.5 split_pcell{si}(gi)+0.5], 'Color', 'white', 'LineWidth', 0.001);
                line([split_pcell{si}(gi)+0.5 split_pcell{si}(gi)+0.5], [0.5 max(split_pcell{si})+0.5], 'Color', 'white', 'LineWidth', 0.001);
            end
            
            
%             if gi == 1
%                     p_start = split_pcell{si}(gi) / 2+0.5;
%             else
%                     p_start = split_pcell{si}(gi-1) +  (split_pcell{si}(gi) - split_pcell{si}(gi-1)) / 2+0.5;
%             end
%             text(p_start, max(split_pcell{si})+ 1.5, label_name{pfm_resort_list(gi)}, 'HorizontalAlignment', 'right', ...
%                 'VerticalAlignment', 'middle', 'Rotation', 45, 'FontSize', 8, 'FontWeight', 'normal');
        end
        
        
        hold off
        print(['/home/jinlong/VDisk2/Jinlong/2024_homologous_parcellation_result/HFIP_result/HFIP_ret/pic_dir/pfm_' sk_list{si} '_' pred_klist{ki} '.eps'],'-depsc','-r600');% save eps with 600dpi

