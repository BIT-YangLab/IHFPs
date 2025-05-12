function eta_to_template_vox = eta2_similarity(X, Y)
% DESCRIPTION:
%   Calculate the similarity value of eta2
%
% USAGE: 
%   X = zeros(100, 1); Y = zeros(100, 1);
%   eta_to_template_vox = eta2_similarity(X, Y)
%
% Inputs:   X,Y: matrix of BOLD signal 
%
% Outputs:  eta_to_template_vox: similarity value between X and Y
%
% X : 80 * 1483
% Y : 17 * 1483
if size(X, 2) ~= size(Y, 2)
    error('dimension not match');
end
m = size(X, 1);
t = size(Y, 1);

eta_to_template_vox = zeros(m, t);
for i = 1:m
    for j = 1:t
        cmap = X(i,:)';
        
        tmap = Y(j,:)';
        Mgrand  = (mean(mean(tmap)) + mean(mean(cmap)))/2;
        Mwithin = (tmap+cmap)/2;
        SSwithin = sum(sum((tmap-Mwithin).*(tmap-Mwithin))) + sum(sum((cmap-Mwithin).*(cmap-Mwithin)));
        SStot    = sum(sum((tmap-Mgrand ).*(tmap-Mgrand ))) + sum(sum((cmap-Mgrand ).*(cmap-Mgrand )));
        eta_to_template_vox(i,j) = 1 - SSwithin/SStot;

    end

end
end