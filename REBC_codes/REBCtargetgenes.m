%targetgenes of the mirs we got from TCGA stuff !
clear; clc;
REBC_full = readtable("gene_REBCsample.csv","ReadVariableNames", true, "ReadRowNames", true);
REBC_tableData = readtable("REBCsample_and_gene_tabledata.csv", "ReadVariableNames", true);
%%
REBC_nozero = REBC_full(5:end, :);
REBC_full_vals = table2array(REBC_nozero);
h = height(REBC_full_vals);
for i = 1:h
    row = h-i+1;
    if all(REBC_full_vals(row,:) == 0)
        REBC_nozero(h-i+1,:) = [];
    end
    disp([string(i) "out of 60664"]);
end
writetable(REBC_nozero, "REBC_test/REBC_data/REBCgene_cleanednozero.csv", ...
    "WriteVariableNames",true, "WriteRowNames", true);

%% 
clear; clc;

relevantMiR = readtable("miR_and_target_genes.xlsx", "ReadVariableNames",true);
REBC_fullgenes = readtable("REBCgene_cleanednozero.csv","ReadVariableNames", true);
relevantMiRgenes = unique(string(table2array(relevantMiR(:,2))));
rebcGenes = string(table2array(REBC_fullgenes(:,1)));

%%
clear geneIndices;
geneIndices = [];
for i = 1:length(relevantMiRgenes)
    gene = relevantMiRgenes(i,1);
    genInd = find(rebcGenes==gene);
    geneIndices = [geneIndices; genInd];
end

geneIndices = sort(geneIndices);

rebcRelevantGenes = REBC_fullgenes(geneIndices, :);

%%
writetable(rebcRelevantGenes, 'REBC_test/REBC_data/REBC_targetgenes_ofTCGAmir.csv','WriteVariableNames',true);
%% no duplicate gene reads
% clc; clear;
% REBCmirTargetGenes = readtable('REBC_targetgenes_ofTCGAmir.csv', 'ReadVariableNames', true);
% REBCgenes = string(table2array(REBCmirTargetGenes(:,1)));
% 
% %find the duplicates in our target genes:
% [uniqueList,~,uniqueNdx] = unique(REBCgenes);
% N = histc(uniqueNdx,1:numel(uniqueList));
% dupNames = uniqueList(N>1);
% dupNdxs = arrayfun(@(x) find(uniqueNdx==x), find(N>1), ...
%     'UniformOutput',false);
% 
% %%
% 
% clear canRemoveInd;
% canRemoveInd = [];
% cantRemoveInd = [];
% for i = 1:height(dupNdxs)
%     dupNdx1 = dupNdxs{i,1};
%     dupNdxvals = table2array(REBCmirTargetGenes(dupNdx1, 2:end));
% %     scatter(dupNdxvals(1,:), dupNdxvals(2,:));
%     x = corrcoef([dupNdxvals(1,:)', dupNdxvals(2,:)']);
%     y = 1;
%     if length(dupNdx1) > 2
%         z = corrcoef([dupNdxvals(1,:)', dupNdxvals(3,:)']);
%         y = z(2,1);
%     end
%     if x(2,1) > 0.98 && x(2,1) < 1.02 && y>0.98 &&y<1.02
%         canRemoveInd = [canRemoveInd; {dupNdx1(2:end,1)}];
%     else
%         cantRemoveInd = [cantRemoveInd; {dupNdx1}];
%     end
% end
% 
% %% what do we do with indices we cant remove :-( rename them haha
% 
% for i = 1:height(cantRemoveInd)
%     dupNdx2 = cantRemoveInd{i,1};
%     geneA = string(REBCmirTargetGenes{dupNdx2(1,1),"gene_name"});
%     geneB = string(REBCmirTargetGenes{dupNdx2(2,1),"gene_name"});
%     REBCmirTargetGenes(dupNdx2(1,1),"gene_name") = {char(strcat(geneA, "_read_v1"))};
%     REBCmirTargetGenes(dupNdx2(2,1),"gene_name") = {char(strcat(geneB, "_read_v2"))};
% end
% 
