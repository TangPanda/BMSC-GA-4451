clear; clc;
%%
%table about data:
tableData = necessarySampleData('gdc_sample_sheet_rna.csv');
%%
geneTPM = organizeTCGA('gdc_sample_sheet_rna.csv');
%%
writetable(geneTPM, 'THCA_raw_gene_TPM.csv', WriteVariableNames=true);
%%
writetable(tableData, 'THCA_gene_about.csv');
