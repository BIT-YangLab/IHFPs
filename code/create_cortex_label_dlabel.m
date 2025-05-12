function create_cortex_label_dlabel(label, outfile, cmp)
% DESCRIPTION:
%   Visualization of atlas with label
%
% USAGE: 
%   cmp = zeros(17, 3);
%   create_cortex_label_dlabel(label, outfile, cmp)
%
% Inputs:   label: label of atlas 
%
%           outfile: result file
%
%           cmp: color map
%
% Outputs:  
%
% Written by Jinlong Li, Guoyuan Yang
    addpath('./cifti-matlab-master');
    template_dlabel = [ getenv('RESOURCE_DIR') '/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_netassignments_LR.dlabel.nii'];
    colormap_file_tmp = [ getenv('RESOURCE_DIR') '/Temp_networks_colormap.txt'];
    
    fid = fopen(colormap_file_tmp, 'w');
    for row_i = 1:size(cmp, 1)
        fwrite(fid, sprintf('%d\n', row_i));
        fwrite(fid, sprintf('%d %d %d %d %d\n', row_i, ceil(cmp(row_i, 1) * 255), ceil(cmp(row_i, 2)*255), ceil(cmp(row_i, 3)*255), 255));
    end
    fclose(fid);
    
    
    template_file_nii = [ getenv('RESOURCE_DIR') '/MNI_immidiate_file.dlabel.nii'];
    template_d = cifti_read(template_dlabel);
    load( [ getenv('RESOURCE_DIR') '/fs_LR_32k_no_medial_wall_index.mat']);  %no_medial_wall_index

    if length(label) == 96854
        label(setdiff(1:64984, no_medial_wall_index)) = [];
    end
    if length(label) ~= 91282
        error('dimension = 91282 or 96854');
    end
    
    new_label = label ;
    

    template_d.cdata = single(new_label);
    cifti_write(template_d, template_file_nii);
    system(['wb_command -cifti-label-import ' template_file_nii ' ' colormap_file_tmp ' ' outfile]);
    system(['rm ' template_file_nii]);
    system(['rm ' colormap_file_tmp]);
end