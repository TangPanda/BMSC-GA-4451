clear; clc;

data_table = readtable("REBC_mirtargetgenes_tumor_TDS.csv");

geneNames = string(data_table.Properties.VariableNames);
thyroid_function_genes = string({'DIO1', 'DIO2', 'DUOX1', 'DUOX2', ...
    'FOXE1', 'GLIS3', 'NKX2_1', 'PAX8', 'SLC26A4', 'SLC5A5', 'SLC5A8', ...
    'TG', 'THRA', 'THRB', 'TPO', 'TSHR'});
TDSgenes = intersect(geneNames,thyroid_function_genes);

data_table_noTDS = removevars(data_table, [TDSgenes(1,1), TDSgenes(1,2)]);
w = width(data_table_noTDS);
obs = data_table_noTDS(:,1);
vars = data_table_noTDS.Properties.VariableNames(2:w-1);
x = table2array(data_table_noTDS(:, 2:w-1));
TDS = table2array(data_table_noTDS(:,w));
%% find best alpha - it is .1

bestalpha = findBestAlpha(x,TDS,vars,50,.1,10);

figure(1)
bar(bestalpha(:,1), bestalpha(:,2));
title('Alpha vs MSE');
xlabel('Alpha');
ylabel('MSE');
xlim([.05, 1.05]);
ylim([0, .057]);
%%
bcoeftable50 = loopRegressCoefs(x,TDS,vars,50,0.1,10);
bcoef50_vars = bcoeftable50.Properties.VariableNames;

%%
clear stableCountCol;
stableCountCol = zeros(height(bcoeftable50), 1);
for i = 1:height(bcoeftable50)
    notzeros = length(find(bcoeftable50{i,:}));
    stability = notzeros/length(bcoef50_vars);
    stableCountCol(i,1) = stability;
end

bcoeftable50.stablecount = stableCountCol;
bcoeftable50.MeanCoef = mean(bcoeftable50{:,bcoef50_vars}, 2);
bcoeftable50.STDcoef = std(bcoeftable50{:,bcoef50_vars}, 0,2);
%% write a table with onlye the genes and expressions and coefs we care abt
h = length(vars);
bcoefhm50 = find(bcoeftable50.MeanCoef(1:h) ~= 0);
bcoefhmImptgenes50 = cell2table(vars(bcoefhm50)');
bcoefhmImptgenes50.meanCoef = bcoeftable50.MeanCoef(bcoefhm50);
bcoefhmImptgenes50.stdCoef = bcoeftable50.STDcoef(bcoefhm50);
bcoefhmImptgenes50.stability = bcoeftable50.stablecount(bcoefhm50);
%%
importantgenes_absweight = maxk(abs(bcoeftable50.MeanCoef(1:h)), 25);

clear importantGenes;
for i = 1:25
    absweight = importantgenes_absweight(i,1);
    index = find(abs(bcoeftable50.MeanCoef) == absweight);
    importantGenes(i,1) = {index};
    importantGenes(i,2) = {vars(index)};
    importantGenes(i,3) = {bcoeftable50.MeanCoef(index)};
    importantGenes(i,4) = {bcoeftable50.STDcoef(index)};
    importantGenes(i,5) = {bcoeftable50.stablecount(index)};
end
%%
figure(3)
bplot = bar(cell2mat(importantGenes(:,3)));
xtik = 1:25;
set(gca, 'XTickLabel', string(importantGenes(:,2)), 'XTick', 1:25);
title('Coefficients of Top 25 Genes in Prediciting TDS - 1SE (\alpha = .1, 50x repeat)');
ylabel('Coefficient of Variable');
xlabel('Genes');
hold on 
er = errorbar(xtik, cell2mat(importantGenes(:,3)),cell2mat(importantGenes(:,4)));  
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off
%%
figure(5)
for i = 1:6
    subplot(2,3,i);
    lassoGene = importantGenes{i,2};
    scatter(TDS, table2array(data_table_noTDS(:, lassoGene)));

    set(gca,'fontsize',16);
    title(['TDS vs' lassoGene 'Expression']);
    xlim([-4.1,4]);
    xlabel('TDS');
    ylabel([lassoGene 'expression']);
    legend('off');
end

%% plot predicted vs actual using lasso regression 10x cv w/ best alpha
clear new_yfit indexes_choose;
new_yfit = ones(length(TDS), 1);
indexes_choose = [1:390];

for element_id = 1:10 % 10fold cv
    clear new_X new_Y indexlist new_index;
    a = randperm(length(indexes_choose),39);
    indexlist = indexes_choose(a);
    new_index = [1:390];
    new_index(indexlist) = [];
    new_Y  = TDS(new_index);
    new_X = x(new_index,:);
    indexes_choose(a)= [];
    [new_B,new_FitInfo] = lasso(new_X,new_Y,'CV',10, 'Alpha', .3);
    new_lam = new_FitInfo.Index1SE;
    new_coefs = new_B(:,new_lam);
    new_intercept = new_FitInfo.Intercept(new_lam);
    new_yfit(indexlist,1) = (x(indexlist,:)*new_coefs) + new_intercept;
end
hold off;
figure(7) %visualize with 10fold cv
hold on;
for row = 1:length(TDS)
    plot(TDS(row),new_yfit(row),'linestyle','none','marker','o','MarkerSize',7, MarkerFaceColor='#D95319');
end
 plot([-5 5],[-5 5],'linestyle','-','linewidth',1,'color','k');
set(gca,'fontsize',20);
title('10-fold cross-validation 1SE (\alpha = .1)');
xlabel('Measured TDS');
ylabel('Predicted TDS');
hold off;
%%
writetable(bcoefhmImptgenes50, 'REBC_test/REBC_data/REBC_50x1seLassoGenes_noTDSgenes_bestalpha.csv');