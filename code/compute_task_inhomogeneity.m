function ret = compute_task_inhomogeneity(label, contrast_map)
% DESCRIPTION:
%   Calculate the functional inhomogeneity 
%
% USAGE: 
%   ret = compute_task_inhomogeneity(label, contrast_map)
%
% Inputs:   label: label of atlas 
%
%           contrast_map: task contrast map
%
%
% Outputs:  ret: value of task inhomogeneity
%
% Written by Jinlong Li, Guoyuan Yang
    label_list = unique(label);
    label_list(label_list < 1) = [];
    label_list(isnan(label_list)) = [];
    ret_std = [];
%     index1 = find(~isnan(contrast_map));
%     s_series = contrast_map(index1);
%     s_series = bsxfun(@minus, s_series, mean(s_series, 1));
%     s_series = bsxfun(@times, s_series, 1./sqrt(sum(s_series.^2, 1)));
%     contrast_map(index1) = s_series;
%     contrast_map(index1) = FisherTransform(contrast_map(index1));
    for l_i = 1:length(label_list)
        index = find(label == label_list(l_i));
        ret_std(end+1) = std(contrast_map(index)) * length(index);
    end
    
    ret = sum(ret_std) / nnz(label);%length(label);
end