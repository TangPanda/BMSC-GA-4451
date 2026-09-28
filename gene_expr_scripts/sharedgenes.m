%% this section only needs to be run once
clear; clc;

%get the genes expressed used in TCGA database
TCGA_full = readtable('obsSamples_varsTCGAgenes.csv', 'ReadVariableNames',true);
TCGA_genes = string(table2array(TCGA_full(:, 1)));

%get the genes expressed in CCLE 
CCLE_full = readtable('thyroid_gene_expression.csv','ReadVariableNames',true);
CCLE_genes = CCLE_full.Properties.VariableNames;

CCLE_genes = string(CCLE_genes);
CCLE_genes_nonum = eraseBetween(CCLE_genes, "_", "_");
CCLE_genes_nonum = (erase(CCLE_genes_nonum, "__"))';

%% this is to make the sheets 
clear both_genes
both_genes = intersect(CCLE_genes_nonum, TCGA_genes);
writematrix(both_genes, 'same_genes.csv')

CCLE_full.Properties.VariableNames = CCLE_genes_nonum;
writetable(CCLE_full,'thyroid_gene_expression_nonum.csv','WriteVariableNames',true);


%% now write tables that only have the shared genes
clear; clc;

shared_genes = readtable('same_genes.csv');
TCGA_full = readtable('cleanednozero_TCGA_genes.csv', 'ReadVariableNames',true);
CCLE_full = readtable('thyroid_gene_expression_nonum.csv','ReadVariableNames',true);


%%
TCGA_genes = string(table2array(TCGA_full(:,1)));
clear THCAindex;
clear CCLE_sharedgenes;
clear h;

h = height(shared_genes);

CCLE_sharedgenes = CCLE_full(:,1);
col = 2;
for i = 1:h
    gene = string(shared_genes{i,1});
    CCLE_sharedgenes(:,col) = CCLE_full(:,char(gene));
    CCLE_sharedgenes.Properties.VariableNames(col) = gene;
    ind = find(strcmp(TCGA_genes,gene)==1);
    THCAindex{i,1} = ind;
    col = col+1;
end


%%
% writetable(CCLE_sharedgenes, 'CCLEthyroid_sharedgenes.csv', 'WriteVariableNames',true);

