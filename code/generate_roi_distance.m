function dis_map = generate_roi_distance(cortical_roi)
% DESCRIPTION:
%   Calculate the Euclidean distance between two roi
%
% USAGE: 
%   cortical_roi = zeros(64984, 1);
%   dis_map = generate_roi_distance(cortical_roi)
%
% Inputs:   cortical_roi: atlas
%
% Outputs:  dis_map:  Euclidean distance between ROIs
%


% ROI center vertices on cortical sheet
cortical_list = unique(cortical_roi);
cortical_list(cortical_list < 1) = [];
cortical_list(isnan(cortical_list)) = [];
cort_parcel_num = length(cortical_list);

lh_cortical_list = unique(cortical_roi(1:32492));
lh_cortical_list(lh_cortical_list < 1) = [];
lh_cortical_list(isnan(lh_cortical_list)) = [];
cort_parcel_num_lh = length(lh_cortical_list);
cort_parcel_num_rh = cort_parcel_num - cort_parcel_num_lh;

% % lh_roi = gifti('./atlas/MBNA_124_32k_L.label.gii');
% lh_roi = gifti('./atlas/CHARM2_SARM2_in_Yerkes19.L.label.gii');
% lh_roi = lh_roi.cdata;
% % rh_roi = gifti('./atlas/MBNA_124_32k_R.label.gii');
% rh_roi = gifti('./atlas/CHARM2_SARM2_in_Yerkes19.R.label.gii');
% rh_roi = rh_roi.cdata;
% % rh_roi(rh_roi~=0) = rh_roi(rh_roi~=0) + 124;
% cortical_roi = [lh_roi; rh_roi];  64984 * 1

roi_center_ind = zeros(cort_parcel_num, 1);

inflated_surf_lh = gifti('/home/jinlong/VDisk2/Jinlong/resources//Conte69.L.very_inflated.32k_fs_LR.surf.gii');
inflated_surf_rh = gifti('/home/jinlong/VDisk2/Jinlong/resources//Conte69.R.very_inflated.32k_fs_LR.surf.gii');
vertices_pos = [inflated_surf_lh.vertices; inflated_surf_rh.vertices];




for i = 1: cort_parcel_num
    roi_vertices_ind = find(cortical_roi==cortical_list(i));
    roi_vertices_pos = vertices_pos(roi_vertices_ind, :);

    dist_x = repmat(roi_vertices_pos(:, 1), 1, size(roi_vertices_pos, 1)) - repmat(roi_vertices_pos(:, 1), 1, size(roi_vertices_pos, 1))';
    dist_y = repmat(roi_vertices_pos(:, 2), 1, size(roi_vertices_pos, 1)) - repmat(roi_vertices_pos(:, 2), 1, size(roi_vertices_pos, 1))';
    dist_z = repmat(roi_vertices_pos(:, 3), 1, size(roi_vertices_pos, 1)) - repmat(roi_vertices_pos(:, 3), 1, size(roi_vertices_pos, 1))';

    dist_W = sqrt(dist_x.^2 + dist_y.^2 + dist_z.^2);

    [~, min_sumdist_pos] = min(sum(dist_W, 1));
    
    roi_center_ind(i) = roi_vertices_ind(min_sumdist_pos);
end

roi_center_ind_lh = roi_center_ind(1:cort_parcel_num_lh);
roi_center_ind_rh = roi_center_ind(cort_parcel_num_lh+1:end);

midthickness_surf_lh_path = '/home/jinlong/VDisk2/Jinlong/resources//Conte69.L.midthickness.32k_fs_LR.surf.gii';
midthickness_surf_rh_path = '/home/jinlong/VDisk2/Jinlong/resources//Conte69.R.midthickness.32k_fs_LR.surf.gii';

result_geo_dist_lh_path = '/home/jinlong/VDisk2/Jinlong/resources//geo_dist_lh.func.gii';
result_geo_dist_rh_path = '/home/jinlong/VDisk2/Jinlong/resources//geo_dist_rh.func.gii';


% iterate
dist_lh = zeros(cort_parcel_num_lh);
dist_rh = zeros(cort_parcel_num_rh);
dist_lh_rh = zeros(cort_parcel_num_lh, cort_parcel_num_rh);
dist_rh_lh = zeros(cort_parcel_num_rh, cort_parcel_num_lh);

for i = 1: cort_parcel_num_lh
    pos_lh = roi_center_ind_lh(i) - 1;
    system(['wb_command -surface-geodesic-distance ', midthickness_surf_lh_path, ' ', int2str(pos_lh), ' ', result_geo_dist_lh_path]);
    dist_lh_i = gifti(result_geo_dist_lh_path);
    dist_lh(i, :) = dist_lh_i.cdata(roi_center_ind_lh);
    dist_lh_rh(i, :) = dist_lh_i.cdata(roi_center_ind_rh - 32492);

    
end

for i = 1 : cort_parcel_num_rh
    pos_rh = roi_center_ind_rh(i) - 1 - 32492;
    system(['wb_command -surface-geodesic-distance ', midthickness_surf_rh_path, ' ', int2str(pos_rh), ' ', result_geo_dist_rh_path]);
    dist_rh_i = gifti(result_geo_dist_rh_path);
    dist_rh(i, :) = dist_rh_i.cdata(roi_center_ind_rh - 32492);
    dist_rh_lh(i, :) = dist_rh_i.cdata(roi_center_ind_lh);

end

dist_lh = (dist_lh + dist_lh') ./ 2;
dist_rh = (dist_rh + dist_rh') ./ 2;

dist_lhrh = (dist_lh_rh + dist_rh_lh') ./ 2;

dis_map = zeros(cort_parcel_num);
dis_map(1:cort_parcel_num_lh, 1:cort_parcel_num_lh) = dist_lh;
dis_map(1:cort_parcel_num_rh, 1:cort_parcel_num_rh) = dist_rh;
dis_map(1:cort_parcel_num_lh, cort_parcel_num_lh+1:cort_parcel_num) = dist_lhrh;
dis_map(cort_parcel_num_lh+1:cort_parcel_num, 1:cort_parcel_num_lh) = dist_lhrh';


end