function no_medial_wall_index =  find_cortical_label_index(mesh)
% DESCRIPTION:
%   remove the index of medial wall
%
% USAGE: 
%   mesh = 'fs_LR_32k';
%   no_medial_wall_index =  find_cortical_label_index(mesh)
%
% Inputs:   mesh: 'fs_LR_32k' or 'fsaverage'
%
% Outputs:  no_medial_wall_index: the index of cortex after removing the medial wall
%
    if ~exist('mesh', 'var')
        mesh = 'fs_LR_32k';
    end
    addpath([ getenv('CBIG_CODE_DIR') '/utilities/matlab/fslr_matlab' ]);
    addpath([ getenv('CBIG_CODE_DIR') '/external_packages/SD/SDv1.5.1-svn593/BasicTools/' ]);

    if contains(mesh, 'fs_LR')
        lh_targ_mesh = CBIG_read_fslr_surface('lh', mesh,'inflated','medialwall.annot');
        rh_targ_mesh = CBIG_read_fslr_surface('rh', mesh,'inflated','medialwall.annot');
        no_medial_wall_index = [find(lh_targ_mesh.MARS_label == 2); find(rh_targ_mesh.MARS_label == 2) + length(lh_targ_mesh.MARS_label)];
    elseif contains(mesh, 'fsaverage') 
        lh_targ_mesh = CBIG_ReadNCAvgMesh('lh', mesh, 'inflated', 'cortex');
        rh_targ_mesh = CBIG_ReadNCAvgMesh('rh', mesh, 'inflated', 'cortex');
        
        no_medial_wall_index = [find(lh_targ_mesh.MARS_label == 2) find(rh_targ_mesh.MARS_label == 2) + length(lh_targ_mesh.MARS_label)]';
    end
end
