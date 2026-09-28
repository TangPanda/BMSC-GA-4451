function renamed_table = renameDuplicateGenes(dupNdxs, full_table)

n = length(dupNdxs);
for i = 1:n
    c = length(dupNdxs{i,1});
    indexarray = dupNdxs{i,1};
    for h = 1:c
        index = indexarray(h,1);
        gene_id = string(table2array(full_table(index, "gene_name")));
        full_table(index, "gene_name") =  {char(strcat(gene_id, "_read_v", string(h)))};
    end
end

renamed_table = full_table;

end


