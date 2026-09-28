clear; clc;

data_matrix_TDS = readtable("mirTargetGenes_TDS.csv", 'ReadRowNames',true, 'ReadVariableNames',true);

obs_TDS = data_matrix_TDS.Properties.RowNames;
vars_TDS = data_matrix_TDS.Properties.VariableNames;
w = length(vars_TDS);
data_matrix_TDS_vals = data_matrix_TDS{:,1:w-1};
%this is TDS scores:
data_matrix_TDS_scores = data_matrix_TDS{:,w};
%PCA!
[coefficientsTDS,scoresTDS,latentTDS,tsquaredTDS,explainedTDS,muTDS] = pca(data_matrix_TDS_vals);
%% now PCA plot!
figure(2);
plot(0:length(explainedTDS),[0 cumsum(explainedTDS)'],'linestyle','-','linewidth',2,'color','k','marker','o','markersize',8,'markerfacecolor','k','markeredgecolor','k');
set(gca,'XLim',[0 30],'YLim',[0 100],'Xtick',0:length(explainedTDS),'Xticklabel',0:length(explainedTDS),'Box','off','fontsize',20);
title('Cumulative variance explained by PCA');
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
figure(3);
for counter = 1:3
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
%%
figure(4)
scatter3(scoresTDS(:,1),scoresTDS(:,2),scoresTDS(:,3), [], markerColors, 'filled',  MarkerEdgeColor='black', SizeData=60, LineWidth=1)
axis equal
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
zlabel('3rd Principal Component')
