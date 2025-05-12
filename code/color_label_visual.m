function color_label_visual(source, out_dir, file_stem, flg, netw)
% DESCRIPTION:
%   Visualization of parcellations
%
% USAGE: 
%   color_label_visual(labels, './', '_test', 64984, 7);
%
% Inputs:   source : labels of parcellation map
%
%           out_dir: directory of result
%        
%           file_stem: used to create the file name
%
%           flg: numbero of vertex
%
%           netw: number of parcels
% Outputs:  
%
% Written by Jinlong Li, Guoyuan Yang
    resources_dir = getenv("RESOURCE_DIR");
    load( [resources_dir '/fs_LR_32k_no_medial_wall_index.mat']);  %no_medial_wall_index
    create_label_script='./Convert_to_cifti_label.py';

    %% color file
    if netw == 7
        color_L_file = [ resources_dir '/Lh_7networks_fs32_new.label.gii'];
        color_R_file = [ resources_dir '/Rh_7networks_fs32_new.label.gii'];
    elseif netw == 14
        color_L_file = [ resources_dir  '/Lh_gordon333_13net.label.gii'];
        color_R_file = [ resources_dir  '/Rh_gordon333_13net.label.gii'];
    elseif netw == 17
        color_L_file = [ resources_dir  '/Lh_17networks_fs32_new_color.label.gii'];
        color_R_file = [ resources_dir  '/Rh_17networks_fs32_new_color.label.gii'];
    else
        color_L_file = [ resources_dir  '/Lh_Gordon_networks_fs32_new_color_fix.label.gii'];
        color_R_file = [ resources_dir  '/Rh_Gordon_networks_fs32_new_color_fix.label.gii'];
    end
 
    if flg == 59412
        new_label = zeros(64984, 1);
        new_label(no_medial_wall_index) = source;
    elseif flg == 64984
        new_label = source;
    end

    lh_labels=new_label(1:32492);
    rh_labels=new_label(32493:64984);
    save([out_dir '/temp_label.mat'], 'lh_labels', 'rh_labels');

    L_out_file = [out_dir '/Gordon_infomap_fill_color_L' file_stem '.label.gii'];
    R_out_file = [out_dir  '/Gordon_infomap_fill_color_R' file_stem '.label.gii'];

    system(['python  ' create_label_script '   '  out_dir '/temp_label.mat  ' color_L_file '  ' color_R_file '  '  L_out_file '  '  R_out_file]);

    
    system(['wb_command -set-structure ' L_out_file ' CORTEX_LEFT']);
    system(['wb_command -set-structure ' R_out_file ' CORTEX_RIGHT']);

end