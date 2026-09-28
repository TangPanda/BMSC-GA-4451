clear; clc;

mir_target_genes_sheet = readtable('hsa_MTI.xlsx', 'ReadVariableNames',true);

miR_tGene = mir_target_genes_sheet(:, ["miRNA" "TargetGene"]);
%%
miR_list = string(table2array(miR_tGene(:,1)));

topMIR = readtable("miRNA_side/miR_data/miR_VIPs_top.xlsx", 'ReadVariableNames',true);

clear indices;
indices = [];
l = height(topMIR);
for i = 1:l
    miR = replace(string(topMIR{i,1}), "_","-");
    miRind = find(contains(miR_list,miR));
    indices = [indices; miRind];
end

indices = sort(indices);

relevantMiR = miR_tGene(indices, :);
% %%
% writetable(relevantMiR,'gene_expr_side/gene_expr_data/miR_and_target_genes.xlsx','WriteVariableNames',true);
%%
TCGA_full = readtable('cleanednozero_TCGA_genes.csv', 'ReadVariableNames',true);
%%
relevantMiRgenes = unique(string(table2array(relevantMiR(:,2))));
tcgaGenes = string(table2array(TCGA_full(:,1)));
%%
clear geneIndices;
geneIndices = [];
for i = 1:length(relevantMiRgenes)
    gene = relevantMiRgenes(i,1);
    genInd = find(tcgaGenes==gene);
    geneIndices = [geneIndices; genInd];
end

geneIndices = sort(geneIndices);

tcgaRelevantGenes = TCGA_full(geneIndices, :);

%%
writetable(tcgaRelevantGenes, 'gene_expr_data/TCGAtarget_genes2.csv','WriteVariableNames',true);
%%
clc; clear;
TCGAmirTargetGenes = readtable('TCGAtarget_genes2.csv', 'ReadVariableNames', true);
TCGAgenes = string(table2array(TCGAmirTargetGenes(:,1)));

%find the duplicates in our target genes:
[uniqueList,~,uniqueNdx] = unique(TCGAgenes);
N = histc(uniqueNdx,1:numel(uniqueList));
dupNames = uniqueList(N>1);
dupNdxs = arrayfun(@(x) find(uniqueNdx==x), find(N>1), ...
    'UniformOutput',false);

%%

clear canRemoveInd;
canRemoveInd = [];
cantRemoveInd = [];
for i = 1:height(dupNdxs)
    dupNdx1 = dupNdxs{i,1};
    dupNdxvals = table2array(TCGAmirTargetGenes(dupNdx1, 2:end));
%     scatter(dupNdxvals(1,:), dupNdxvals(2,:));
    x = corrcoef([dupNdxvals(1,:)', dupNdxvals(2,:)']);
    y = 1;
    if length(dupNdx1) > 2
        z = corrcoef([dupNdxvals(1,:)', dupNdxvals(3,:)']);
        y = z(2,1);
    end
    if x(2,1) > 0.98 && x(2,1) < 1.02 && y>0.98 &&y<1.02
        canRemoveInd = [canRemoveInd; {dupNdx1(2:end,1)}];
    else
        cantRemoveInd = [cantRemoveInd; {dupNdx1}];
    end
end

%% what do we do with indices we cant remove :-( rename them haha

for i = 1:height(cantRemoveInd)
    dupNdx2 = cantRemoveInd{i,1};
    geneA = string(TCGAmirTargetGenes{dupNdx2(1,1),"gene_name"});
    geneB = string(TCGAmirTargetGenes{dupNdx2(2,1),"gene_name"});
    TCGAmirTargetGenes(dupNdx2(1,1),"gene_name") = {char(strcat(geneA, "_read_v1"))};
    TCGAmirTargetGenes(dupNdx2(2,1),"gene_name") = {char(strcat(geneB, "_read_v2"))};
end

%% now lets remove the indices that are duplicates
canRemoveInd = sort(cell2mat(canRemoveInd));
h = length(canRemoveInd);
TCGAmirTargetGenesTable = TCGAmirTargetGenes;
for i = 1:h
    indToBeRemoved = canRemoveInd(h-i+1,1);
    TCGAmirTargetGenesTable(indToBeRemoved,:) = [];
end
%% 


geneNames = string(table2array(TCGAmirTargetGenesTable(:,1)));
uniq = unique(geneNames);
TCGAmirTargetGenesTable = transposeTable(TCGAmirTargetGenesTable);
%%
writetable(TCGAmirTargetGenesTable,'gene_expr_data/mirtargetgenes_datamatrix2.csv', 'WriteVariableNames',true, 'WriteRowNames',true);