function homo_pro_prior =  find_homologous_area(label_list)
%% find_homogeneity_area(indlabel_list)
% using overlap region to find the most homogeneity area
% indlabel_list -> matrix  label_cnt * sub_num : label of each individual as the column, label:64984*1

%% prepare
if isempty(label_list)
    error('label_list is empty');
end

%% create prior prob
homo_pro_prior = zeros(64984,400);
for pa = 1:400
    for i = 1:size(label_list,2)
        homo_pro_prior(:, pa) = homo_pro_prior(:, pa) + (label_list(:,i) == pa);
    end
end

homo_pro_prior = homo_pro_prior ./ size(label_list,2);

end