clear; clc;

coefs50xgenes = readtable('50x1seLassoGenes_noTDSgenes.csv', 'ReadVariableNames',true);
coefs50xgenes.Properties.VariableNames(1) = {'Gene'};
genelist_1se = string(coefs50xgenes.Gene);
mir_genetarget = readtable('miR_and_target_genes.xlsx');
%% just 1se miR first
allmirgenetar = string(mir_genetarget.TargetGene);
clear mirtable;
for i = 1:length(genelist_1se)
    gene = genelist_1se(i,1);
    idx = find(string(mir_genetarget.TargetGene)==gene);
    mirtable(1:length(idx),i) = mir_genetarget.miRNA(idx);
end
mirtable = cell2table(mirtable);
mirtable.Properties.VariableNames = genelist_1se;

%%
clear mirtable_otherway;
topmirs = string(table2array(readtable('miRNA_side/miR_data/miR_VIPs_top 2.xlsx', 'ReadVariableNames', false)));
allmir = replace(string(mir_genetarget.miRNA), "-","_");
allmir = replace(allmir, "hsa_", "");
for i = 1:length(topmirs)
    mir = topmirs(1,i);
    idxm = find(allmir == mir);
    targenes = intersect(unique(string(mir_genetarget.TargetGene(idxm))), genelist_1se);
    mirtable_otherway(1:length(targenes),i) = targenes;
end
mirtable_otherway = array2table(mirtable_otherway);
mirtable_otherway.Properties.VariableNames = topmirs;
%% write the tables
% writetable(mirtable_otherway, 'gene_expr_side/gene_expr_data/50xmir_imptgenes.csv');
% writetable(mirtable, 'gene_expr_side/gene_expr_data/imptgenes_mir2.csv');

%% I want to plot mirna expression vs gene expression for the top 5 genes!
clear; clc;
mir_genestable = readtable('50xmir_imptgenes.csv', 'ReadRowNames', false);
genes_tdsTable = readtable('TCGAgenes_nozscore_TDS.csv');
mir_expressionTable = readtable('miRNA_side/miR_data/TCGAmiRNAlog2PTC.xlsx');
%%
genes_samples = string(genes_tdsTable.Row);
mir_Table_inorder = rearrange(mir_expressionTable, genes_samples);
mir_samples = string(mir_Table_inorder.obs);
removesample = setdiff(genes_samples, mir_samples);
for i = 1:length(removesample)
    ind = find(string(genes_tdsTable.Row) == removesample(length(removesample)-i+1));
    genes_tdsTable(ind,:) = [];
end
%% plot mirs vs gene w tds as color >.<
clear datamatrix markerColors;
datamatrix = genes_tdsTable.TDS;
markerColors = addTDSColor(datamatrix,150);
%%
for w = 1:width(mir_genestable)
    figure(w) %1st mirna
    mirname = string(mir_genestable.Properties.VariableNames(w));
    mirgenelist = mir_genestable(:,w);
    h = height(mir_genestable);
    for i = 1:h
        if string(mir_genestable.(mirname)(h-i+1)) == ""
            mirgenelist(h-i+1,:) = [];
        end
    end
    for j = 1:height(mirgenelist)
        subplot(ceil(height(mirgenelist)/2),2,j);
        gene = char(mirgenelist{j,1});
        geneexpr = genes_tdsTable.(gene);
        mir_expr = mir_Table_inorder.(mirname);
        a = scatter(geneexpr,mir_expr, [], markerColors, 'filled', ...
            MarkerEdgeColor='black', SizeData=60, LineWidth=1);
        set(gca,'fontsize',10);
        title([gene ' expression VS ' newline char(replace(mirname, "_", "-")) ' expression']);
        xlabel([gene ' expression' newline '(TPM, log2, z-scored)']);
        ylabel([char(replace(mirname, "_", "-")) 'expression'...
            newline '(RPM, log2, z-scored)']);
        legend('off');
    end
    sgtitle([char(replace(mirname, "_", "-")) ' versus Target Genes with TDS as color'])
    colormap(redblue(150))
    cbh = colorbar(); 
    caxis([0,1])
    set(cbh, 'YTick', [0, 1], ...
        'YTickLabel',{'low TDS', 'High TDS'})
end

%%
genes_mirtable = readtable('imptgenes_mir.csv', 'ReadRowNames', false);
topgenes = ["IGFBPL1"; "KIT"; "SLC1A1"; "LIPG"; "SCUBE3"; "DOCK5"; "TNS4"; ...
    "DPYSL5"; "DPT"; "TPPP"];

for v = 1:length(topgenes)
    figure(v+9) %1st gene
    genename = topgenes(v,1);
    genemirlist = genes_mirtable.(genename);
    h = height(genes_mirtable);
    for i = 1:h
        if string(genes_mirtable.(genename)(h-i+1)) == ''
            genemirlist(h-i+1,:) = [];
        end
    end
    genemirlist = unique(string(genemirlist));
    geneexpr = genes_tdsTable.(genename);
    for j = 1:length(genemirlist)
        subplot(ceil(length(genemirlist)/2),2,j);
        mir = replace(string(genemirlist{j,1}), "-", "_");
        mir = replace(mir,"hsa_","");
        mir_expr = mir_Table_inorder.(mir);
        scatter(geneexpr,mir_expr, [], markerColors, 'filled', ...
            MarkerEdgeColor='black', SizeData=60, LineWidth=1);
        set(gca,'fontsize',10);
        title([char(genename) ' vs. ' char(replace(mir, "_", "-")) ' with TDS as color']);
        xlabel([char(genename) ' expression' "(TPM, log2 normalized)"]);
        ylabel([char(replace(mir, "_", "-")) 'expression' "(RPM, log2 normalized)"]);
        legend('off');
    end
    colormap(redblue(150))
    cbh = colorbar(); 
    caxis([0,1])
    set(cbh, 'YTick', linspace(0,1,9), ...
        'YTickLabel',string(linspace(-4,4,9)))
    cbh.Title.String = "TDS";
end
%%
figure(1)
scatter(genes_tdsTable.TDS, genes_tdsTable.DPYSL5);



