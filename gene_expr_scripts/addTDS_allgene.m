clear; clc;

%will this even work... heres to hoping

TDS_scores = readtable("miRNA_side/miR_data/TCGA_name_and_TDS_PTC.xlsx", 'ReadVariableNames',true);
%i am addin tds to notz-scored, only log2, but u can change the file name
%depending on what u wanna add the tds scores to
data_table_noTDS = readtable('obsSamples_varsTCGAgenes_log2.csv', 'ReadRowNames',true);
%%
colTDS = array2table(zeros(571,1) - 2000);
data_table_noTDS(:, 57364) = colTDS;

TDSnames_inorder = replace(string(table2array(TDS_scores(:,1))), "-", "_");
data_table_names_inorder = string(data_table_noTDS.Properties.RowNames);

dataTableTDS = data_table_noTDS;

for i = 1:length(TDSnames_inorder)
    TDSname = TDSnames_inorder(i,1);
    if isnan(find(strcmp(data_table_names_inorder, TDSname))) == 0
        dataTableTDS(TDSname, 57364) = TDS_scores(i,2);
    end
end

%%
h = height(dataTableTDS);
dataTableTDS_fin = dataTableTDS;
for i = 1:h
    if dataTableTDS{h-i+1,57364} == -2000
        dataTableTDS_fin(h-i+1,:) = [];
    elseif isnan(dataTableTDS{h-i+1,57364})
        dataTableTDS_fin(h-i+1,:) = [];
    end
end

%%
dataTableTDS_fin.Properties.VariableNames(57364) = {'TDS'};
%%
writetable(dataTableTDS_fin, 'gene_expr_data/TCGAgenes_nozscore_TDS.csv', 'WriteRowNames',true,'WriteVariableNames',true);
