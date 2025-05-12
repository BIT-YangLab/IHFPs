function [homo_ci, std_ci] = compute_parcellation_mean_fc(labels, input_filename)
% DESCRIPTION:
%   Calculate the functional measures of individual atlas
%
% USAGE: 
%   labels = zeros(64984, 1);
%   input_filename = {'/path/to/cifti/file'};
%   [homo_ci, std_ci] = compute_parcellation_mean_fc(labels, input_filename)
%   global_mean_all = homo_ci(:, 1);
%   global_mean_within = homo_ci(:, 2);
%   global_mean_between = homo_ci(:, 3);
%
% Inputs:   labels: parcellation maps
%
%           input_filename: cifti file
%
%
% Outputs:  homo_ci: Nx3, N parcels columns: global mean fc, within-system fc, between-system fc
%
%           std_ci: Nx3, N parcels columns: std of global mean fc, within-system fc, between-system fc
% Written by Jinlong Li, Guoyuan Yang
    if ~exist(input_filename)
        error('file not exist');
    end

    labels(isnan(labels)) = 0;
    

    cifti_s = ft_read_cifti(input_filename);
    dt = cifti_s.dtseries(1:64984, :);
    corr1 = paircorr_mod(dt');
    corr1(corr1<0) = 0;


    label_t_list = unique(labels);
    label_t_list(label_t_list < 1) = [];
    label_t_list(isnan(label_t_list)) = [];
    homo_ci = zeros(length(label_t_list), 3);
    std_ci = zeros(length(label_t_list), 3);

    for li1 = 1:length(label_t_list)
        idx1 = find(labels == label_t_list(li1));
        within_fc = corr1(idx1, idx1); 
        within_fc = within_fc(~eye(size(within_fc, 1)));
        homo_ci(li1, 1) = nanmean(within_fc(:));
        std_ci(li1, 1) = nanstd(within_fc(:));

        between_fc = corr1(idx1, setdiff(1:size(corr1, 2), idx1));
        homo_ci(li1, 2) = nanmean(between_fc(:));
        std_ci(li1, 2) = nanstd(between_fc(:));

        all_fc = [within_fc(:); between_fc(:)];
        homo_ci(li1, 3) = nanmean(all_fc(:));
        std_ci(li1, 3) = nanstd(all_fc(:));

    end
end