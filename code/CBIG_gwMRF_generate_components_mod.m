function [lh_ci,lh_sizes]=CBIG_gwMRF_generate_components_mod(label,roi_nei)

    % This function computes connected components for a labeled vector.
    %
    % - Input
    %   - lh_avg_mesh = left hemisphere mesh
    %   - rh_avg_mesh = right hemisphere mesh
    %   - lh_labels = label vector for left hemisphere for corresponding mesh
    %   - rh_labels = label vector for right hemisphere for corresponding mesh
    %
    % Ouput
    %   - lh_ci = connected components for the left hemisphere
    %   - rh_ci = connected components for the right hemisphere
    %   - lh_sizes = size of the corressponding components for the left hemisphere
    %   - rh_sizes = size of the corressponding components for the right hemisphere
    %
    % Example
    %   - [lh_ci,rh_ci,lh_sizes,rh_sizes]=CBIG_gwMRF_generate_components(lh_avg_mesh,rh_avg_mesh,lh_labels,rh_labels)
    %Written by Alexander Schaefer and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    % original code: if node a and node b have the same label, set d = 1
    tic
    
    nk = size(roi_nei, 1);
    G = zeros(size(roi_nei));
    idx = find(roi_nei);
    for j_idx = 1:length(idx)
        idx_i = idx(j_idx);
        rk = floor((idx_i-1) / nk) + 1;
        ck = mod(idx_i-1, nk) + 1; 
        G(rk, ck) = label(rk) == label(ck);
    end

    g = sparse(G);
    toc
    [lh_ci, lh_sizes] = components(g);
    
       
    
end
