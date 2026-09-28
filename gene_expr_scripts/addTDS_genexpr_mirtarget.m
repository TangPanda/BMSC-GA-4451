clear; clc;

TDS_scores = readtable("miRNA_side/miR_data/TCGA_name_and_TDS_PTC.xlsx", 'ReadVariableNames',true);

data_table_noTDS = readtable('zscored_mirtargetgenes2.csv', 'ReadRowNames',true);
%%
colTDS = array2table(zeros(571,1) - 2000);
data_table_noTDS(:, 2070) = colTDS;

TDSnames_inorder = replace(string(table2array(TDS_scores(:,1))), "-", "_");
data_table_names_inorder = string(data_table_noTDS.Properties.RowNames);

dataTableTDS = data_table_noTDS;

for i = 1:length(TDSnames_inorder)
    TDSname = TDSnames_inorder(i,1);
    if isnan(find(strcmp(data_table_names_inorder, TDSname))) == 0
        dataTableTDS(TDSname, 2070) = TDS_scores(i,2);
    end
end

%%
h = height(dataTableTDS);
dataTableTDS_fin = dataTableTDS;
for i = 1:h
    if dataTableTDS{h-i+1,2070} == -2000
        dataTableTDS_fin(h-i+1,:) = [];
    elseif isnan(dataTableTDS{h-i+1,2070})
        dataTableTDS_fin(h-i+1,:) = [];
    end
end

%%
dataTableTDS_fin.Properties.VariableNames(2070) = {'TDS'};
%%
writetable(dataTableTDS_fin, 'mirTargetGenes_TDS_2.csv', 'WriteRowNames',true,'WriteVariableNames',true);
