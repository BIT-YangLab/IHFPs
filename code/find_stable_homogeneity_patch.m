function stable_net = find_stable_homogeneity_patch(new_label_list)
%% find_stable_homogeneity_patch(label_list)
% find the stable patch found in every sub
% new_label_list -> matrix 64984 * subnum : label of all individuals

    patch_num = max(max(new_label_list));
    label_list = zeros(patch_num, 1);
    sub_cnt = 0;

    sub_n = size(new_label_list,2);
    for i = 1:sub_n
        g_l = new_label_list(:, i);
        g_l_mmap = unique(g_l);
        g_l_mmap(g_l_mmap < 1) = [];
        label_list(g_l_mmap) = label_list(g_l_mmap) + 1;
        sub_cnt = sub_cnt + 1;
    end

    stable_net = find(label_list == sub_cnt);
end