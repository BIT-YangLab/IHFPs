function [dice_list, label_list] = compute_dice(label1, label2)
% DESCRIPTION:
%   Calculate the DICE's coefficient between two atlases
%
% USAGE: 
%   [dice_list, label_list] = compute_dice(label1, label2)
%
% Inputs:   label1, labels: two parcellation maps
%
%
% Outputs:  dice_list: list of DICE's coefficient
%
%           label_list: list of labels

    label_list = union(unique(label1), unique(label2));
    label_list(label_list < 1) = [];
    label_list(isnan(label_list)) = [];

    dice_list = zeros(length(label_list), 1);
    for pi = 1:length(label_list)
        index1 = label1 == label_list(pi);
        index2 = label2 == label_list(pi);
        index = index1 + index2;
        dice_list(pi) = nnz(index == 2) / nnz(index);
    end
end