function eta_to_template_vox = template_matching_eta2(template_profile, participant_profile)

% load group parcellation and network 400 * N, 400 1st is group parcel,
% N 2st is reference group parcel(self)

% individual parcellation, 400 * N , 400 is individual parcel, N is group
% ref parcel, this matrix is connectome profile

% using eta2 similarity, 400*N x N*400, -> 400 * 400, get the ind-grp,
% assign network label of each parcel -> 400 * 1

if size(template_profile, 2) ~= size(participant_profile, 2)
    error('dimension not match');
end

eta_to_template_vox = zeros(size(participant_profile, 1), size(template_profile, 1));

for i = 1:size(template_profile, 1)
    for j = 1:size(participant_profile, 1)
        goodvox = (~isnan(template_profile(i,:)) & ~isnan(participant_profile(j, :)));
        cmap = template_profile(i,goodvox);
        %tmap = cifti_template_mat_full(j,goodvox)';
        tmap = participant_profile(j, goodvox);
        Mgrand  = (mean(mean(tmap)) + mean(mean(cmap)))/2;
        Mwithin = (tmap+cmap)/2;
        SSwithin = sum(sum((tmap-Mwithin).*(tmap-Mwithin))) + sum(sum((cmap-Mwithin).*(cmap-Mwithin)));
        SStot    = sum(sum((tmap-Mgrand ).*(tmap-Mgrand ))) + sum(sum((cmap-Mgrand ).*(cmap-Mgrand )));
        eta_to_template_vox(j,i) = 1 - SSwithin/SStot;

        
    end
end


end

