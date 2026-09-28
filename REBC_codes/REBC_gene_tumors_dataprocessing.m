clear; clc;

%cluster genes :-)
data_matrix = readtable('REBC_targetgenes_ofTCGAmir.csv');
data_matrix_right = transposeTable(data_matrix);
%%
data_matrix_log2 = log2tableROWNAME(data_matrix_right);
%%
zScored_data_matrix = zScoreTable(data_matrix_log2,1);
%%
writetable(zScored_data_matrix, 'REBC_test/REBC_data/REBC_zscored_mirtargetgenes_allsamp.csv', 'WriteVariableNames',true);
%% this part is to actually cluster our data!
clear; clc;
data_matrix = readtable('REBC_zscored_mirtargetgenes_allsamp.csv');
observations = data_matrix{:,1};
vars = data_matrix.Properties.VariableNames([2:end]);
data_matrix_vals = data_matrix{:,2:end};
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cluster1 = clustergram(data_matrix_vals,'Cluster',3,'Symmetric',true,'Colormap',redblue(150),...
    'DisplayRange',3,'RowPDist','correlation','ColumnPDist','euclidean','Linkage','complete',...
    'RowLabels',observations,'ColumnLabels',vars); 
%% just get tumor data!
clear; clc;
data_matrix = readtable('REBC_targetgenes_ofTCGAmir.csv');
data_matrix_right = transposeTable(data_matrix);
% %% write
% writetable(data_matrix_right, "REBC_test/REBC_data/REBCsample_targetgenes.csv", "WriteRowNames",true, "WriteVariableNames",true);

tableinfo = readtable("REBCsample_and_gene_tabledata.csv", "ReadVariableNames", true);
%%
datatype = unique(tableinfo.Var6);
idxtumor = find(matches(string(tableinfo.Var6), "Primary Tumor"));
tumorsamples = tableinfo.Var5(idxtumor);
normalsamples = tableinfo.Var5(setdiff(1:end, idxtumor));

clear data_matrix_right;
data_matrix_right = readtable("REBCsample_targetgenes.csv", "ReadVariableNames", true);


tumor_normal = string([tumorsamples; normalsamples]);
index_start_normal = length(tumorsamples)+1;
REBC_tumor_normal = rearrange(data_matrix_right, tumor_normal);

clear tumor_gene_samp normal_gene_samp;
tumor_gene_samp = REBC_tumor_normal(1:index_start_normal-1,:);
normal_gene_samp = REBC_tumor_normal(index_start_normal:end, :);
% %%
% writetable(tumor_gene_samp, "REBC_test/REBC_data/REBCtumorsample_targetgenes.csv", "WriteVariableNames",true);
% writetable(normal_gene_samp, "REBC_test/REBC_data/REBCnormalsample_targetgenes.csv", "WriteVariableNames",true);
%%
%tumor_gene_samp = readtable("REBCtumorsample_targetgenes.csv", ... 
%"ReadVariableNames", true);
tumor_log2 = log2table(tumor_gene_samp);
%%
zScored_tumor = zScoreTable(tumor_log2,1);
%%
writetable(zScored_tumor, 'REBC_test/REBC_data/REBC_zscored_mirtargetgenes_tumor.csv', 'WriteVariableNames',true);
%% this part is to actually cluster our data!
data_matrix = readtable('REBC_zscored_mirtargetgenes_tumor.csv');
observations = data_matrix{:,1};
vars = data_matrix.Properties.VariableNames([2:end]);
data_matrix_vals = data_matrix{:,2:end};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cluster1 = clustergram(data_matrix_vals,'Cluster',3,'Symmetric',true,'Colormap',redblue(150),...
    'DisplayRange',3,'RowPDist','correlation','ColumnPDist','euclidean','Linkage','complete',...
    'RowLabels',observations,'ColumnLabels',vars); 
%% rename and add TDS
clear; clc;
data_matrix = readtable('REBC_zscored_mirtargetgenes_tumor.csv');
tableinfo = readtable("REBCsample_and_gene_tabledata.csv", "ReadVariableNames", true);

wantedNames = string([tableinfo.Var5, tableinfo.Var4]);

renamed_tumorgenes = rename(wantedNames, data_matrix);
% %% write sheet of renamed
% writetable(renamed_tumorgenes, 'REBC_test/REBC_data/renamedREBC_zscored_mirtargetgenes_tumor.csv', 'WriteVariableNames',true);

sample_TDS = readtable("REBCsample_and_TDS.xlsx","ReadVariableNames", true);
renamed_tumor_TDS = addTDS(renamed_tumorgenes, sample_TDS);

%% write sheet
writetable(renamed_tumor_TDS, 'REBC_test/REBC_data/REBC_mirtargetgenes_tumor_TDS.csv', 'WriteVariableNames',true);
%%

