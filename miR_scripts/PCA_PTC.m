clear; clc;

%be in the folder miR_scripts

%% use the miRNA and TDS data to see TDS in color as our '3rd dimension'
data_matrix_TDS = readtable('02_THCA_miR_withTDS_noNormal.xlsx', 'ReadRowNames', true);
%establish observations, variables, and width of our data matrix
obs_TDS = data_matrix_TDS.Properties.RowNames;
vars_TDS = data_matrix_TDS.Properties.VariableNames;
w = length(vars_TDS);
%this is miRNA expression values:
data_matrix_TDS_vals = data_matrix_TDS{:,1:w-1};
%this is TDS scores:
data_matrix_TDS_scores = data_matrix_TDS{:,w};
%PCA!
[coefficientsTDS,scoresTDS,latentTDS,tsquaredTDS,explainedTDS,muTDS] = pca(data_matrix_TDS_vals);



%% now PCA plot!
figure(1);
plot(0:length(explainedTDS),[0 cumsum(explainedTDS)'],'linestyle','-','linewidth',2,'color','k','marker','o','markersize',8,'markerfacecolor','k','markeredgecolor','k');
set(gca,'XLim',[0 30],'YLim',[0 100],'Xtick',0:length(explainedTDS),'Xticklabel',0:length(explainedTDS),'Box','off','fontsize',20);
title('Cumulative variance explained by PCA', 'FontSize', 18);
xlabel('Principal components (PCs)');
ylabel('% Variance explained');

% mapcaplot(data_matrix_TDS_vals);
%% to assign colors to my TDS scores
cmap = colormap(redblue(150));
yourValues = data_matrix_TDS_scores;
v = rescale(yourValues, 1, 150); % Nifty trick!
numValues = length(yourValues);
markerColors = zeros(numValues, 3);
% Now assign marker colors according to the value of the data.
for k = 1 : numValues
    row = round(v(k));
    markerColors(k, :) = cmap(row, :);
end
%%
figure(2)
    scatter(scoresTDS(:,1),scoresTDS(:,2), [], markerColors, 'filled',  MarkerEdgeColor='black', SizeData=60, LineWidth=1)
    set(gca,'fontsize',20);
    title('PCA score plot');
    xlabel('PC1');
    ylabel(['PC2']);
    legend('off');
    colormap(redblue(150))
    cbh = colorbar; 
    % set color range
    caxis([0,1])
    % set ticks 
    set(cbh, 'YTick', [0, 1], ...
        'YTickLabel',{'low TDS', 'High TDS'})

%% scatter plot along PCs 
figure(3);
for counter = 1:4
    subplot(2,2,counter);
    hold on;
    scatter(scoresTDS(:,1),scoresTDS(:,counter+1), [], markerColors, 'filled',  MarkerEdgeColor='black', SizeData=60, LineWidth=1)

    set(gca,'fontsize',20);
    title('PCA score plot');
    xlabel('PC1');
    ylabel(['PC' num2str(counter+1)]);
    legend('off');
    colormap(redblue(150))
    cbh = colorbar; 
    % set color range
    caxis([0,1])
    % set ticks 
    set(cbh, 'YTick', [0, 1], ...
        'YTickLabel',{'low TDS', 'High TDS'})
end

%% 3d PCA plot
figure(4)
scatter3(scoresTDS(:,1),scoresTDS(:,2),scoresTDS(:,3), [], markerColors, 'filled',  MarkerEdgeColor='black', SizeData=60, LineWidth=1)
axis equal
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
zlabel('3rd Principal Component')
%% PCA Coefficient plot (PC1 versus PC3 bar plot)
vars_TDS = replace(string(vars_TDS),"_", "-");
figure(5);
bar(coefficientsTDS(1:92,[1 2 3]));
set(gca,'XTick',1:92,'XTickLabel',vars_TDS,'YLim',[-0.3 0.3]);
title('PCA loadings plot');
xlabel('Signals');
ylabel('PCA loadings');
xtickangle(45);
legend (['PC1'; 'PC2'; 'PC3']);
%% top 5 variables most important in each PC

%how many most important vars do you want?
n_impt = 5;
%how many PCs do you want to look at?
PCnum = 3;

topLoadTable = topNimportant(coefficientsTDS,PCnum,n_impt);

clear topLoadingmiRNA
for i = 1:PCnum
    indexColumn = (i-1)*2+1;
    indices = topLoadTable(:,indexColumn);
    l = height(indices);
    for j = 1:l
        indexmiRNA = table2array(indices(j,1));
        topLoadingmiRNA{j, indexColumn} = vars_TDS(1,indexmiRNA);
    end
    topLoadingmiRNA(:,i*2) = table2cell(topLoadTable(:, i*2));
end
topLoadingmiRNA = cell2table(topLoadingmiRNA);
topLoadingmiRNA.Properties.VariableNames = ["top 5 miRNAs in PC1", "PC1 Loading Values", "top 5 miRNAs in PC2", "PC2 Loading Values", "top 5 miRNAs in PC3", "PC3 Loading Values"];

figure(9)
for i = 1:PCnum
    subplot(1,3,i)
    loadingvals = table2array(topLoadingmiRNA(:,i*2));
    bar(loadingvals);
    set(gca,'XTick',1:n_impt,'XTickLabel',replace(table2array(topLoadingmiRNA(1:n_impt, (i-1)*2+1)), "_", "-"),'YLim',[0 0.3]);
    title(['Top '  num2str(n_impt) ' PC' num2str(i) ' loadings plot (absolute value)']);
    xlabel('miRNA');
    ylabel('|PCA loadings|');
    xtickangle(45);
end

%% get table of top loading miRNAs in #of PCS as figure
fig = figure;
uit = uitable(fig,'Data',table2cell(topLoadingmiRNA));
uit.ColumnName={topLoadingmiRNA.Properties.VariableNames{:,:}};
uit.RowName=[]; %removing default row numbering as in your uitable
