function new_label = create_ref_groupROI_single_group(output_dir, group_label, ind_label_list)
%% create_ref_groupROI(output, group_label, ind_label_list)
% create age_independent group ref-ROI
% output_dir -> str: result ref_ROI file's path
% group_label -> matrix 64984 * 1: age_independent label
% ind_label_list -> cell group_cnt * 1: each group's individual label  

%% prepare
parcel_n_cnt = 400;

if ~exist(output_dir)
    mkdir(output_dir);
end

if length(group_label) ~= 64984
    error('array size wrong: group_label');
end

if size(group_label, 2) == 64984
    group_label = group_label';
end

% if ~isa(ind_label_list, 'cell')
%     error('type of ind_label_list wrong');
% end



%% create group ref ROI




    
    % find homogeneity area 

    homo_pro_prior = find_homologous_area(ind_label_list);

    new_label = zeros(64984,1);
    MuI_threshhold_all_networks = findoverlapthreshold(homo_pro_prior, 1:parcel_n_cnt, 0,  'smooth_then_derivative');
    MuI_threshhold_all_networks(MuI_threshhold_all_networks == 0) = 0.75;
    for i = 1:parcel_n_cnt % parcel number: hard code (modify)
        
        index = find(homo_pro_prior(:, i) >= MuI_threshhold_all_networks(i));
        if isempty(index)
            continue;
        end
        label_match_list = unique(group_label(index));
        label_match_list(label_match_list < 1) = [];

        match_list = zeros(length(label_match_list), 2);
        for j = 1:length(label_match_list)
            index1 = find(group_label == label_match_list(j));
            index2 = intersect(index1, index);
            match_list(j, 1) = length(index2) ;
            match_list(j, 2) = length(index1);
    
        end
    
    
        [~, m_i] = maxk(match_list(:, 1), 2);
        target_p = label_match_list(m_i(1));
        if length(m_i) > 1 && match_list(m_i(1), 1) / match_list(m_i(2), 1) < 2 && match_list(m_i(1), 1) / match_list(m_i(1), 2) < match_list(m_i(2), 1) / match_list(m_i(2), 2)
            target_p = label_match_list(m_i(2));
    
        end
        index3 = find(group_label == target_p);
        new_label(intersect(index3, index)) = target_p;
    end

        new_label_l = new_label;

    



parcel_list = unique(new_label_l);
parcel_list(parcel_list < 1) = [];

for i = 1:length(parcel_list)
    pa = parcel_list(i);
    index = find(new_label_l == pa);
    if length(index) < 6
        new_label_l(new_label_l == pa) = 0;
    end
end

new_label = new_label_l;
save([output_dir '/group_ref_roi.mat'], 'new_label');

end