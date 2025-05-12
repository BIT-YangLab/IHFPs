function info_sub = read_hcd_info()

a = importdata('/home/jinlong/VDisk2/datasets/HCD/Package_1188937_hcpdrestingRecommand/fmriresults01.txt');

index = [];

b = strsplit(a{1}, '"');
index = [];
ig_index = [6 12 15 16 24 26 31 ];
for k = 1:length(b)
    if nnz(isstrprop(b{k}, 'alphanum'))
        index(end+1) = k;
    end
end
head_info = b(index);
head_info(ig_index) = [];



for i = 3:length(a)
    b = strsplit(a{i}, '"');
    index = [];
    for k = 1:length(b)
        if nnz(isstrprop(b{k}, 'alphanum'))
            index(end+1) = k;
        end
    end
    b = b(index);
    for k  =1:length(b)
        if i == 3
            eval(['info_sub.' head_info{k} ' =   b(k) ;' ]);
        else
            eval(['info_sub.' head_info{k} '(end+1) =  b(k) ;' ]);
        end
    end
end


end