function hex_colors = rgb2hex(rgb_matrix)
% DESCRIPTION:
%   Transform RGB to hex code 
%
% USAGE: 
%   rgb_matrix = rand(17, 3);
%   hex_colors = rgb2hex(rgb_matrix)
%
% Inputs:   rgb_matrix: matrix of RGB, Nx3
%
% Outputs:  hex_colors: hex code, N*1
    hex_colors = cell(size(rgb_matrix, 1), 1);
    for i = 1:size(rgb_matrix, 1)
        hex_colors{i} = sprintf('#%02X%02X%02X', rgb_matrix(i, 1), rgb_matrix(i, 2), rgb_matrix(i, 3));
    end
end

