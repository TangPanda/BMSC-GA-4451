%this script is to compare patient survival/tumor phenotype vs the important genes!

clc; clear;

%first read in important genes from our lasso regression model (alpha = .3)
importantGenes = readtable('50x1seLassoGenes_noTDSgenes_bestalpha.csv');
%and we want to read our gene expression dataset
geneExpression = readtable('zscored_mirtargetgenes2.csv');
%finally read the patient data !
patientData = readtable('clinical.tsv', 'FileType', 'text');
%% Get data about top genes and the patient phenotype/classification

% decide how many genes you care about - maybe we will start with 5 
[~, topGenes_idx] = maxk(abs(importantGenes.meanCoef), 6);
topGenes = importantGenes(topGenes_idx, :);

%Now we need to get only the gene expressions of the topGenes!
topGeneXpr = geneExpression(:,["Row"; string(topGenes.Var1)]);

%also lets match the samples from gene expression to patientData
patientData_noXtra = patientData([1:2:height(patientData)], :);


%% now we want to create a table with sample ID, gene expr, phenotype

%table structure:
% | sample ID| gene expr| ... |phenotype info| ... 


%first we match the samples of gene expression and patient data by case id

sampleID = replace(string(table2array(patientData_noXtra(:,"case_submitter_id"))), "-", "_");
pattern = "_" + digitsPattern(2) + lettersPattern(1);

sampleGeneID = replace(string(table2array(topGeneXpr(:,"Row"))), pattern, "");

[bothSampl, idPatientData, idGeneData] = intersect(sampleID, sampleGeneID);
%% make table of data
%what patient info do we want? - make a table of info in order of intersect
%case_submitter_id, age_at_index, days_to_death, vital_status, 
% ajcc_pathologic_m, ajcc_pathologic_n,
%ajcc_pathologic_stage, ajcc_pathologic_t, ajcc_staging_system_edition,
%morphology, primary_diagnosis
vars_of_interest = ["case_submitter_id", "age_at_index", "days_to_death", "vital_status", ...
    "ajcc_pathologic_m", "ajcc_pathologic_n", "ajcc_pathologic_stage", "ajcc_pathologic_t", ...
    "ajcc_staging_system_edition", "morphology", "primary_diagnosis"];

%now we make a new table!
clear genePhenotypeTbl;
genePhenotypeTbl = array2table(bothSampl(:,1));

for i = 1:length(bothSampl)
    patientIDx = idPatientData(i,1);
    genePhenotypeTbl(i, 2:length(vars_of_interest)+1) = patientData_noXtra(patientIDx, vars_of_interest);
    geneIDx = idGeneData(i,1);
    genePhenotypeTbl(i, (length(vars_of_interest)+2):length(vars_of_interest)+width(topGeneXpr)) = topGeneXpr(geneIDx, 2:end);
end

genePhenotypeTbl.Properties.VariableNames = ["Sample ID", vars_of_interest, string(topGeneXpr.Properties.VariableNames(2:end))];
% %% write it into a sheet
% writetable(genePhenotypeTbl, "gene_expr_data/geneXP_phenotype.csv", "WriteVariableNames",true);
%% make a nested table with categorized gene and phenotype counts
%need to make a table with counts for each phenotype 
genes = topGenes.Var1;
phenotypes = ["ajcc_pathologic_n", "ajcc_pathologic_stage"];
genePhenotypeCounts = classifyNest(genePhenotypeTbl, genes, phenotypes);
%% Graph gene vs phenotype 
%we can start with IGFBPL1

for figgy = 1:width(genePhenotypeCounts)
    figure(figgy*2-1)
    gene = genePhenotypeCounts.Properties.VariableNames{figgy};
    for i = 1:height(genePhenotypeCounts)
        subplot(1,2,i)
        %stacked bar plot with raw numbers of pathologic stage in high or low
        bar([1,2], genePhenotypeCounts.(gene){i,1}{:,:}', 'stacked');
        legend(genePhenotypeCounts.(gene){i,1}.Properties.RowNames);
        set(gca, 'xticklabel', {['Low ' gene], ['High ' gene]});
        ylabel('Raw counts (samples)')
        title(replace(string(genePhenotypeCounts.Properties.RowNames{i}), "_","-"));
    end
    sgtitle(['patient phenotype with ' gene ' high versus low expression - raw count']);

%stacked bar plot with percentage of high or low categorized 

    figure(figgy*2)
    for j = 1:height(genePhenotypeCounts)
        subplot(1,2,j)
        sumh = sum(genePhenotypeCounts.(gene){j,1}{:,1});
        suml = sum(genePhenotypeCounts.(gene){j,1}{:,2});
        percentages = [genePhenotypeCounts.(gene){j,1}{:,1}./sumh, genePhenotypeCounts.(gene){j,1}{:,2}./suml].*100;
        bar([1,2], percentages', 'stacked');
        legend(genePhenotypeCounts.(gene){j,1}.Properties.RowNames);
        set(gca, 'xticklabel', {['Low ' gene], ['High ' gene]}, 'YLim', [0, 100]);
        ylabel('Percentage');
        title(replace(string(genePhenotypeCounts.Properties.RowNames{j}), "_","-"));
    end
    sgtitle(['patient phenotype with ' gene ' high versus low expression - percentage']);

end
%% gene expression per phenotype type! percentage!
for k = 1:height(genePhenotypeCounts)
    figure(k)
    phenotype = genePhenotypeCounts.Properties.RowNames{k};
    for g = 1:width(genePhenotypeCounts)
        gene = genePhenotypeCounts.Properties.VariableNames{g};
        subplot(2,3,g)
        suml = sum(genePhenotypeCounts.(gene){k,1}{:,1});
        sumh = sum(genePhenotypeCounts.(gene){k,1}{:,2});
        percentages = [genePhenotypeCounts.(gene){k,1}{:,1}./suml, genePhenotypeCounts.(gene){k,1}{:,2}./sumh].*100;
        bar([1,2], percentages', 'stacked');
        lgd = legend(genePhenotypeCounts.(gene){k,1}.Properties.RowNames);
        lgd.FontSize = 14;
        set(gca, 'xticklabel', {['Low ' gene], ['High ' gene]}, 'YLim', [0, 100]);
        ylabel('Percentage');
        title(gene);
    
    end
    pheno = replace(string(phenotype), "_","-");
    sgtitle([pheno ' Tumor Diagnosis and High versus Low Expression of Predictive Genes (%)']);

end

%% gene expr of those who died

deadTbl = genePhenotypeTbl(find(matches(genePhenotypeTbl.vital_status, 'Dead')), :);
for i = 1:5
    figure(i)
    gene = topGenes.Var1{i};
    histogram(deadTbl.(gene),[-3 0 3]);
    set(gca, 'XTick', [-1.5 1.5] , 'xticklabel', {['Low ' gene], ['High ' gene]});
    ylabel('Number of Dead Patients');
    title(['Division of dead patients between high and low ' gene 'expression']);
end
