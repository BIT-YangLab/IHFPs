function new_label = trans_to_label(ref_label, id_list, value_list)
    new_label = zeros(size(ref_label));
    for ki = 1:size(ref_label, 2)
        for pi = 1:length(id_list)
            idx = find(ref_label == id_list(pi));
            if isempty(idx)
                continue;
            end
            new_label(idx, ki) = value_list(pi);
        end
    end