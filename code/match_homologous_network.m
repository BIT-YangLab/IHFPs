function ind_new_label =  match_homologous_network(ref_label, ind_label_orig, ind_cifti, num_session)
%% match_homologous_network(ref_label, ind_label_orig, ind_cifti)
% match individual patch to group ref_label
% ref_label -> (64984,1): group atlas contain the most homologous ROI
% ind_label_orig -> (64984, 1) : individual parcellation derived from
% Kong2022
% ind_cifti -> cell{num_session, 1} : individual cifti file path
% num_session -> int : number of session

%% check input parameters
if length(ref_label) ~= length(ind_label_orig)
    error('array size cannot match: ref_label && ind_label_orig');
end

if ~isa(ind_cifti, 'cell')
    error('type of ref_label wrong');
end

if length(ind_cifti) ~= num_session
    error('array size cannot match: ind_cifti && num_session');
end


%% main procedure
patch_list = unique(ind_label_orig);
patch_list(patch_list < 1) = [];
ind_new_label = zeros(64984,1);

% read individual cifti file to create FC profile
fc_profile = cell(num_session,1);
for se_i = 1:num_session
    cifti_file_path = ind_cifti{se_i, 1};
    cifti_s = ft_read_cifti(cifti_file_path);
    fc_profile{se_i, 1} = cifti_s.dtseries(1:64984,:);
end

% match procedure
for p_i = 1:length(patch_list)
    
    pat_id = patch_list(p_i);
%     if pat_id == 173
%             pat_id;
%     end
    index = find(ind_label_orig == pat_id);
    g_a_t = ref_label(index);
    g_tmp_list = unique(g_a_t);
    g_tmp_list(g_tmp_list < 1) = [];
    max_v = mode(g_a_t(g_a_t > 0));
    index_2 = find(g_a_t == max_v);
    n_2 = length(index_2);
    if n_2 < 4
        continue;
    end

    % find kernel ref network
    most_con_list = [];
    for j = 1:length(g_tmp_list)
        index_1 = find(g_a_t == g_tmp_list(j));
        n_1 = length(index_1);
        n_3 = length(find(ref_label == g_tmp_list(j)));
%             if n_2 / n_1 < 2 || (n_1 > 10 && n_1/n_3 > 0.5)
%                 most_con_list(end+1) = g_tmp_list(j);
%             end
        if n_1 / n_3 > 0.3
            most_con_list(end+1) = g_tmp_list(j);
        end
    end

    % three issues: one-one, mul-one(same as one-one), one-mul (ind patch - ref ROI)
    cnt_n = length(most_con_list);
    if cnt_n == 1  % one patch match one kernel ROI 
        ind_new_label(index) = most_con_list(1);
        continue;
    elseif cnt_n == 0
        % if one-mul, but multi-ROI not in kernel area, add the marginal,
        % ROI to split the patch 

        % Only one ROI but not bigger than 0.3, so just match
        if length(g_tmp_list) == 1 
            ind_new_label(index) = g_tmp_list(1);
            continue;
        end

        for j = 1:length(g_tmp_list)
            index_1 = find(g_a_t == g_tmp_list(j));
        
            n_1 = length(index_1);
            if n_1 > 4
                most_con_list(end+1) = g_tmp_list(j);
            end
        end
    end

    % split patch to multi ref ROI, use Pearson correlation to asign the
    % label
    cnt_n = length(most_con_list);
    if cnt_n > 1
        % compute FC 
        index = find(ind_label_orig == pat_id);
        for  j = 1:length(fc_profile)
            profile = fc_profile{j,1};
            if j == 1
                fc = paircorr_mod(profile(index, :)');
            else
                fc = fc + paircorr_mod(profile(index)');
            end
        end
        fc = fc ./ length(fc_profile);

        % assign corresponding area label
        fc_new = zeros(size(fc,1), cnt_n);
        for j = 1:cnt_n
            index_1 = find(ref_label == most_con_list(j));
            index_2 = intersect(index, index_1);
            ind_new_label(index_2) = most_con_list(j);

            index_1 = find(g_a_t == most_con_list(j));
            fc_new(:, j) = mean(fc(:, index_1), 2);
        end

        % assign the lasting vertices 
        for j = 1:length(index)
            if ind_new_label(index(j)) ~= 0
                continue
            end
            [~, m_i] = max(fc_new(j, :));
            ind_new_label(index(j)) = most_con_list(m_i);
        end
    end
end



end