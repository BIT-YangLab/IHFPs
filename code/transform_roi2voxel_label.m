function voxel_label = transform_roi2voxel_label(roi_list, roi_label, voxel_list)
    voxel_label = zeros(size(voxel_list, 1), size(roi_label, 2));
    for ri = 1:length(roi_list)
        index = find(voxel_list == roi_list(ri));
        voxel_label(index, :) = repmat(roi_label(ri, :), length(index), 1);
    end
end