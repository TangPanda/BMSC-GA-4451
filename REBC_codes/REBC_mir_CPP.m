%this script is useless but if u want u can run it 
%clear; clc;
% % 
% 
% mirTumor_TDS = readtable("REBC_test/REBC_data/DEmirREBC_tumor_TDS.csv", "ReadVariableNames",true);
% mirNormal_TDS = readtable("REBC_test/REBC_data/DEmirREBC_norm_TDS.csv", "ReadVariableNames",true);
% 
% mirTumor_noTDS = mirTumor_TDS(:, 1:width(mirTumor_TDS)-1);
% mirNormal_noTDS = mirNormal_TDS(:, 1:width(mirNormal_TDS)-1);
% 
% %% cluster tumor mirs
% observations_tumors = mirTumor_noTDS{:,1};
% vars_tumors = mirTumor_noTDS.Properties.VariableNames([2:end]);
% % vars = vars(1,2:end);
% data_matrix_vals_tumor = mirTumor_noTDS{:,2:end};
% clustertumor = clustergram(data_matrix_vals_tumor,'Cluster',3,'Symmetric',true,'Colormap',redblue(150),...
%     'DisplayRange',3,'RowPDist','correlation','ColumnPDist','euclidean','Linkage','complete',...
%     'RowLabels',observations_tumors,'ColumnLabels',vars_tumors); 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% cluster normal mirs
% observations_normal = mirNormal_noTDS{:,1};
% vars_normal = mirNormal_noTDS.Properties.VariableNames([2:end]);
% % vars = vars(1,2:end);
% data_matrix_vals_norm = mirNormal_noTDS{:,2:end};
% clusternormal = clustergram(data_matrix_vals_norm,'Cluster',3,'Symmetric',true,'Colormap',redblue(150),...
%     'DisplayRange',3,'RowPDist','correlation','ColumnPDist','euclidean','Linkage','complete',...
%     'RowLabels',observations_normal,'ColumnLabels',vars_normal); 
% 
% %% PCA mirs
% %% use the miRNA and TDS data to see TDS in color as our '3rd dimension'
% data_matrix_TDS = mirTumor_TDS;
% %establish observations, variables, and width of our data matrix
% obs_TDS = data_matrix_TDS{:,1};
% vars_TDS = data_matrix_TDS.Properties.VariableNames(2:end);
% w = width(data_matrix_TDS);
% %this is miRNA expression values:
% data_matrix_TDS_vals = data_matrix_TDS{:,2:w-1};
% %this is TDS scores:
% data_matrix_TDS_scores = data_matrix_TDS{:,w};
% %PCA!
% [coefficientsTDS,scoresTDS,latentTDS,tsquaredTDS,explainedTDS,muTDS] = pca(data_matrix_TDS_vals);
% 
% 
% 
% %% now PCA plot!
% figure(1);
% plot(0:length(explainedTDS),[0 cumsum(explainedTDS)'],'linestyle','-','linewidth',2,'color','k','marker','o','markersize',8,'markerfacecolor','k','markeredgecolor','k');
% set(gca,'XLim',[0 30],'YLim',[0 100],'Xtick',0:length(explainedTDS),'Xticklabel',0:length(explainedTDS),'Box','off','fontsize',20);
% title('Cumulative variance explained by PCA');
% xlabel('Principal components (PCs)');
% ylabel('% Variance explained');
% 
% % mapcaplot(data_matrix_TDS_vals);
% %% to assign colors to my TDS scores
% cmap = colormap(redblue(150));
% yourValues = data_matrix_TDS_scores;
% v = rescale(yourValues, 1, 150); % Nifty trick!
% numValues = length(yourValues);
% markerColors = zeros(numValues, 3);
% % Now assign marker colors according to the value of the data.
% for k = 1 : numValues
%     row = round(v(k));
%     markerColors(k, :) = cmap(row, :);
% end
% %% scatter plot along PCs 
% figure(2);
% for counter = 1:4
%     subplot(2,2,counter);
%     hold on;
%     scatter(scoresTDS(:,1),scoresTDS(:,counter+1), [], markerColors, 'filled',  MarkerEdgeColor='black', SizeData=60, LineWidth=1)
% 
%     set(gca,'fontsize',20);
%     title('PCA score plot');
%     xlabel('PC1');
%     ylabel(['PC' num2str(counter+1)]);
%     legend('off');
%     colormap(redblue(150))
%     cbh = colorbar; 
%     % set color range
%     caxis([0,1])
%     % set ticks 
%     set(cbh, 'YTick', [0, 1], ...
%         'YTickLabel',{'low TDS', 'High TDS'})
% end
% %% PLSR 
% clear X Y;
% X = data_matrix_TDS_vals;
% Y = data_matrix_TDS_scores; %response vars, our TDS
% 
% clear TSS;
% clear PLSR_XLoading PLSR_YLoading PLSR_XScore PLSR_YScore PLSR_yfit;
% clear Rsquare Qsquare;
% 
% TSS = sum((Y-mean(Y)).^2);
% ncomp = min(87,length(Y)-1);
% clear XLoading YLoading XScore YScore BETA PCTVAR MSE stats;
% [XLoading,YLoading,XScore,YScore,BETA,PCTVAR,MSE,stats] = plsregress(X,Y,ncomp,'cv',length(Y));
% 
% % Prediction accuracy (leave-one-out cross validation)
% Qsquare = [0 1-length(Y)*MSE(2,2:end)/TSS];
% % Performance
% Rsquare = [0 cumsum(PCTVAR(2,:))];
% 
% figure(3); % R2 and Q2 evaluation
% plot(0:ncomp,100*Rsquare,'-b','linewidth',2,'marker','o','markersize',10');
% hold on;
% plot(0:ncomp,100*Qsquare,'-r','linewidth',2,'marker','o','markersize',10');
% set(gca,'YLim',[0 100],'Box','off','XTick',0:1:20, 'XLim', [0 20]);
% xlabel('PLS Component');
% ylabel('% Variance Explained/Predicted');
% legend({'R^2' 'Q^2'});
% set(gca,'fontsize',20);
% hold off;
% 
% %% Display model-predicted data versus measured (observed) data using leave-one-out strategy
% 
% % optimum number of PLS component
% optimized_ncomp = 3;
% 
% % Rerun the models using the optimized number of PLS components
% clear new_yfit;
% for element_id = 1:length(Y) % pick an element, remove it and predict it using the rest of the dataset
%     clear new_X new_Y;
%     if element_id == 1
%         new_Y = Y(element_id+1:end);
%         new_X = X(element_id+1:end,:);
%     elseif element_id == length(Y)
%         new_Y = Y(1:element_id-1);
%         new_X = X(1:element_id-1,:);
%     else
%         new_Y = [Y(1:element_id-1); Y(element_id+1:end)];
%         new_X = [X(1:element_id-1,:); X(element_id+1:end,:)];
%     end
%     [new_XLoading,new_YLoading,new_XScore,new_YScore,new_BETA,new_PCTVAR,new_MSE,new_stats] = plsregress(new_X,new_Y,optimized_ncomp);
%     new_yfit(element_id) = [ones(size(X(element_id,:),1),1) X(element_id,:)]*new_BETA;
% end
% hold off;
% figure(4); % Visualize correlation between measured and predicted responses (through leave-one-out cross-validation)
% hold on;
% row = 0;
% for row = 1:length(Y)
%     plot(Y(row),new_yfit(row),'linestyle','none','marker','o','MarkerSize',7, MarkerFaceColor='#D95319');
% end
%  plot([-5 5],[-5 5],'linestyle','-','linewidth',1,'color','k');
% set(gca,'fontsize',20);
% title('Leave-one-out cross-validation');
% xlabel('Measured TDS');
% ylabel('Predicted TDS');
% hold off;
% 
% %% Calculate VIP (Variable Importance in Projection) Scores
% sum1 = zeros(1,92);
% sum2 = 0;
% clear SS Wnorm2;
% for i = 1:optimized_ncomp
%     SS(i) = (YLoading(i)^2)*(XScore(:,i)'*XScore(:,i));
% end
% for i = 1:optimized_ncomp
%     sum2 = sum2 + SS(i);
%     Wnorm2(i) = stats.W(:,i)'*stats.W(:,i);
% end
% 
% clear VIP;
% for counter = 1:92
%     for k = 1:optimized_ncomp
%         sum1(counter) = sum1(counter) + SS(k)*stats.W(counter,k)^2/Wnorm2(k);
%     end
%     VIP(counter) = (92*sum1(counter)/sum2)^0.5;
% end
% 
% % Plot VIP scores
% figure(5);
% bar(VIP, 'FaceColor', '#cfa7fa');
% set(gca,'XTick',1:92,'XTickLabel',vars_TDS,'YLim',[0 3.6],'Fontsize',10);
% title('Variable Importance in Projection');
% xlim=get(gca,'xlim');
% hold on
% plot(xlim,[1 1])
% plot(xlim,[1.5 1.5], 'cyan')
% xlabel('miRNA');
% ylabel('VIP score');
% xtickangle(45);
% hold off;
% 
% %% table of most vip vars w scores over 1.5
% vipOver1 = find(VIP(:,:)>=1);
% varOver1 = vars_TDS(1,vipOver1);
% varOver1vals = VIP(1, vipOver1);
% 
