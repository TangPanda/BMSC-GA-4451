clear; clc;

full_sheet = readtable("cleanednozero_TCGA_genes.csv");

full_sheet = full_sheet(5:end, :);
genes = string(table2array(full_sheet(:,1)));

%find the duplicates in our target genes:
[uniqueList,~,uniqueNdx] = unique(genes);
N = histc(uniqueNdx,1:numel(uniqueList));
dupNames = uniqueList(N>1);
dupNdxs = arrayfun(@(x) find(uniqueNdx==x), find(N>1), ...
    'UniformOutput',false);
%%

renamed_genetable = renameDuplicateGenes(dupNdxs, full_sheet);
renamed_Tgenetable = transposeTable(renamed_genetable);
%%
writetable(renamed_Tgenetable, 'obsSamples_varsTCGAgenes.csv',...
    'WriteRowNames',true,'WriteVariableNames',true);
%%
% clear; clc;

data_table = readtable('obsSamples_varsTCGAgenes.csv');
%%
data_table_log2 = log2table(data_table);
%%
writetable(data_table_log2,'obsSamples_varsTCGAgenes_log2.csv', 'WriteRowNames',true, 'WriteVariableNames',true);
%%
data_table_zscore = zScoreTable(data_table_log2,1);
%%
writetable(data_table_zscore, 'obsSamples_varsTCGAgenes_log2zscore.csv','WriteRowNames',true, 'WriteVariableNames',true);
