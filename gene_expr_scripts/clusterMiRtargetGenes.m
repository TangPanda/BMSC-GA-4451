clear; clc;

%cluster genes :-)
data_matrix = readtable('mirtargetgenes_datamatrix2.csv');

%%
data_matrix_log2 = log2table(data_matrix);
%%
zScored_data_matrix = zScoreTable(data_matrix_log2,1);
%%
writetable(zScored_data_matrix, 'gene_expr_data/zscored_mirtargetgenes2.csv', 'WriteVariableNames',true);
%% this part is to actually cluster our data!
clear; clc;
data_matrix = readtable('zscored_mirtargetgenes2.csv');
observations = data_matrix{:,1};
vars = data_matrix.Properties.VariableNames([2:end]);
data_matrix_vals = data_matrix{:,2:end};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cluster1 = clustergram(data_matrix_vals,'Cluster',3,'Symmetric',true,'Colormap',redblue(150),...
    'DisplayRange',3,'RowPDist','correlation','ColumnPDist','euclidean','Linkage','complete',...
    'RowLabels',observations,'ColumnLabels',vars); 
