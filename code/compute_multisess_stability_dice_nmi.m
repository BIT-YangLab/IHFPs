

load([getenv('OUT_RET_DIR') '/IHFP_ret/IHFP_group_network_large_orig_homo.mat']);
load([getenv('OUT_RET_DIR') '/IHFP_ret/IHFP_parcel_idx.mat']);


% figure 1b
    %% compute dice coef.
    load([getenv('OUT_RET_DIR') '/IHFP_ret/IHFP_multisess_individual_parcellations.mat']);
    
    
    process_n = {'IHFP_label_list', 'Kong2022_label_list'};
    for pi = 1:0
        eval(['new_label_list = ' process_n{pi} ';']);
        
        subn = size(new_label_list{1}, 2);
        intra_l = zeros(400, 591); intra_l(:) = nan;
        inter_l = zeros(400, 591*590/2); inter_l(:) = nan;

        nc_inter = 0;
        for si = 1:subn-1
            disp(['  ssi: ' num2str(si) '/590' ]);
            for sj = si+1:subn
                [v1, l1] = compute_dice(new_label_list{1}(:, si), new_label_list{1}(:, sj));
                nc_inter = nc_inter+1;
                inter_l(l1, nc_inter) = v1;
            end
        end

        for si = 1:subn
                disp(['  ssi: ' num2str(si) '/591' ]);
                [v1, l1] = compute_dice(new_label_list{1}(:, si), new_label_list{2}(:, si));
                intra_l(l1, si) = v1;
        end

        inter_dice = nanmean(inter_l, 2);
        intra_dice = nanmean(intra_l, 2);
        inter_dice(filter_parcel_idx) = 0;
        intra_dice(filter_parcel_idx) = 0;

        save([getenv('OUT_RET_DIR') '/IHFP_ret/stability_dice_' process_n{pi} ], 'inter_l', 'intra_l', 'inter_dice', 'intra_dice');
    end
    
    %% compute NMI
    for pi = 1:2
        eval(['new_label_list = ' process_n{pi} ';']);
        
        s1_l = new_label_list{1};
        s2_l = new_label_list{2};
        subn = size(new_label_list{1}, 2);

        inter = [];
        intra = [];
        for ki = 1:subn-1
            for kj = ki+1:subn
                la = s1_l(:, ki);la(isnan(la)) = 0;
                lb = s1_l(:, kj); lb(isnan(lb)) = 0;
                z = compute_nmi(la, lb);
                inter(end+1) = z;
            end
        end
        for ki = 1:subn
            la = s1_l(:, ki);la(isnan(la)) = 0;
            lb = s2_l(:, ki); lb(isnan(lb)) = 0;
            z = compute_nmi(la, lb);
            intra(end+1) = z;
        end

        figure;
        color1 = [0 0.4470 0.7410];
        color2 = [0.8500 0.3250 0.0980];
        [counts, ~] = histcounts(inter, 0.1:0.0025:1);
        counts = counts ./ length(inter);
        [counts1, ~] = histcounts(intra, 0.1:0.0025:1);
        counts1 = counts1 ./ length(intra);
        b1 = bar(0.101:0.0025:1, counts, 'FaceColor', color1, 'EdgeColor', 'none');
        hold on;
        b2 = bar(0.101:0.0025:1, counts1, 'FaceColor', color2, 'EdgeColor', 'none');
        alpha(b1, 0.6);
        alpha(b2, 0.6);
        xlim([0.85 1])
        % legend('inter', 'intra');
        ax = gca;
        ax.Box = 'off';
        xticklabels([])
        yticklabels([])

         save([getenv('OUT_RET_DIR') '/IHFP_ret/stability_nmi_' process_n{pi} ], 'inter', 'intra');
    end