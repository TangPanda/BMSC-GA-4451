%REBC data file and processing 
%first TDS and sample ID

clear; clc;
%% write data sheet with sample ID and respective TDS scores
sampleinfotable = readtable("raw_data/REBC_THYR/abg2538-data-s1.txt", "ReadVariableNames",true);

TDS_sample_table = sampleinfotable(:, ["REBC_ID", "TDS"]);

writetable(TDS_sample_table, "REBC_test/REBC_data/REBCsample_and_TDS.xlsx", WriteVariableNames=true);
%% write data sheet with sample ID and resp miRNAs 

mirTableData = necessarySampleData('REBC_samplesheet_mir_wnormal.csv');
mirREBC = organizeTCGA('REBC_samplesheet_mir_wnormal.csv');
%%
writetable(mirTableData, "REBC_test/REBC_data/REBCsample_and_mir_tabledata.xlsx", WriteVariableNames=true);
writetable(mirREBC, "REBC_test/REBC_data/mir_REBCsample.csv", WriteVariableNames=true);

%% write data sheet wit sample ID and resp genes

geneTableData = necessarySampleData('REBC_samplesheet_rna_wnormal.csv');
geneREBC = organizeTCGA('REBC_samplesheet_rna_wnormal.csv');
%%
writetable(geneTableData,"REBC_test/REBC_data/REBCsample_and_gene_tabledata.csv", WriteVariableNames=true);
writetable(geneREBC, "REBC_test/REBC_data/gene_REBCsample.csv", WriteVariableNames=true);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MIR SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% separate the samples and the mirs by normal and tumor samples (mir)
clear; clc;
mirTableData = readtable('REBCsample_and_mir_tabledata.xlsx', "ReadVariableNames", true);
mirREbC = readtable('mir_REBCsample.csv', "ReadVariableNames", true);
datatype = unique(mirTableData.Var6);
idxtumor = find(matches(string(mirTableData.Var6), "Primary Tumor"));
% %% get the samples that are tumor and normal
tumormirsamples = mirTableData.Var5(idxtumor);
normalmirsamples = mirTableData.Var5(setdiff(1:end, idxtumor));
% %% write into sheets
% writecell(tumormirsamples,"REBC_test/REBC_data/tumormirsamples.csv");
% writecell(normalmirsamples,"REBC_test/REBC_data/normalmirsamples.csv");


%% tranpose and normalize the mir table for rebc (log2)

t_mirREBC = transposeTable(mirREbC);
log2REBCmir = log2tableROWNAME(t_mirREBC);
%% now split into tumor and non-tumor! 
%make a column of tumor on top, normal on bottom:
tumor_normal = string([tumormirsamples; normalmirsamples]);
index_start_normal = length(tumormirsamples)+1;
REBC_tumor_normal_mir = rearrange(log2REBCmir, tumor_normal);

clear tumor_mirs_l2 normal_mirs_l2;
tumor_mirs_l2 = REBC_tumor_normal_mir(1:index_start_normal-1,:);
normal_mirs_l2 = REBC_tumor_normal_mir(index_start_normal:end, :);
% writetable(tumor_mirs_l2, "REBC_test/REBC_data/mirtumor_log2.csv", "WriteVariableNames",true);
% writetable(normal_mirs_l2, "REBC_test/REBC_data/mirnorm_log2.csv", "WriteVariableNames",true);

%% and zscore write into sheets
zL2_tumormirs = zScoreTable(tumor_mirs_l2, 1);
zL2_normmirs = zScoreTable(normal_mirs_l2, 1);
% %write into sheets
% writetable(zL2_tumormirs, "REBC_test/REBC_data/lessmirtumor_l2zscored.csv", "WriteVariableNames",true);
% writetable(zL2_normmirs, "REBC_test/REBC_data/lessmirnorm_l2zscored.csv", "WriteVariableNames",true);
%% rename the samples and then add TDS
wantedNames = string([mirTableData.Var5, mirTableData.Var4]);

renamed_tumormirs = rename(wantedNames, zL2_tumormirs);
renamed_normmirs = rename(wantedNames, zL2_normmirs);

% %% write renamed into sheets
% writetable(renamed_tumormirs, "REBC_test/REBC_data/less_renamed_mirtumor_l2zscored.csv", "WriteVariableNames",true);
% writetable(renamed_normmirs, "REBC_test/REBC_data/less_renamed_mirnorm_l2zscored.csv", "WriteVariableNames",true);

%%
namesTDSTable = readtable("REBC_test/REBC_data/REBCsample_and_TDS.xlsx");
tumorwithTDS = addTDS(renamed_tumormirs, namesTDSTable);
normwithTDS = addTDS(renamed_normmirs, namesTDSTable);

% %% write into sheets
% writetable(tumorwithTDS,"REBC_test/REBC_data/lessmirREBC_tumor_TDS.csv", WriteVariableNames=true);
% 
% writetable(normwithTDS,"REBC_test/REBC_data/lessmirREBC_norm_TDS.csv", WriteVariableNames=true);
% 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GENE SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ill get back to this haha
clear; clc;
namesTDSTable = readtable("REBC_test/REBC_data/REBCsample_and_TDS.xlsx");
geneREBC = readtable('gene_REBCsample.csv', "ReadVariableNames", true);
geneTableData = readtable('REBCsample_and_gene_tabledata.csv', "ReadVariableNames", true);