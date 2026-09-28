clear; clc;

%add the TDS data to see 
data_matrix_TDS = readtable('02_THCA_miR_withTDS_noNormal.xlsx', 'ReadRowNames', true);
obs_TDS = data_matrix_TDS.Properties.RowNames;
vars_TDS = data_matrix_TDS.Properties.VariableNames;
% vars = vars(1,:);
w = length(vars_TDS);
data_matrix_TDS_vals = data_matrix_TDS{:,1:w-1};
data_matrix_TDS_scores = data_matrix_TDS{:,w};
[coefficientsTDS,scoresTDS,latentTDS,tsquaredTDS,explainedTDS,muTDS] = pca(data_matrix_TDS_vals);

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

%% PLSR
clear X Y;
X = data_matrix_TDS_vals;
Y = data_matrix_TDS_scores; %response vars, our TDS

clear TSS;
clear PLSR_XLoading PLSR_YLoading PLSR_XScore PLSR_YScore PLSR_yfit;
clear Rsquare Qsquare;

TSS = sum((Y-mean(Y)).^2);
ncomp = min(87,length(Y)-1);
clear XLoading YLoading XScore YScore BETA PCTVAR MSE stats;
[XLoading,YLoading,XScore,YScore,BETA,PCTVAR,MSE,stats] = plsregress(X,Y,ncomp,'cv',10);

% Prediction accuracy (leave-one-out cross validation)
Qsquare = [0 1-length(Y)*MSE(2,2:end)/TSS];
% Performance
Rsquare = [0 cumsum(PCTVAR(2,:))];

figure(1); % R2 and Q2 evaluation
plot(0:ncomp,100*Rsquare,'-b','linewidth',2,'marker','o','markersize',10');
hold on;
plot(0:ncomp,100*Qsquare,'-r','linewidth',2,'marker','o','markersize',10');
set(gca,'YLim',[0 100],'Box','off','XTick',0:1:20, 'XLim', [0 20]);
title('PLSR Variance explained and predicted using 10-fold cross validation')
xlabel('PLS Component');
ylabel('% Variance Explained/Predicted');
legend({'R^2' 'Q^2'});
set(gca,'fontsize',20);
hold off;

%% Display model-predicted data versus measured (observed) data using 10fold cv strategy

% optimum number of PLS component
optimized_ncomp = 4;
clear new_yfit indexes_choose;
indexes_choose = [1:height(data_matrix_TDS)];
new_yfit = ones(length(Y), 1);

% Rerun the models using the optimized number of PLS components
clear new_yfit;
for element_id = 1:10 %10fold cv
    clear new_X new_Y indexlist new_index;
    
    if ismember(element_id, [1:mod(height(data_matrix_TDS),10)])
    a = randperm(length(indexes_choose),floor(height(data_matrix_TDS)/10)+1);
    indexlist = indexes_choose(a);
    new_index = [1:height(data_matrix_TDS)];
    new_index(indexlist) = [];
    new_Y  = Y(new_index);
    new_X = X(new_index,:);
    indexes_choose(a)= [];
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
    else
    a = randperm(length(indexes_choose),floor(height(data_matrix_TDS)/10));
    indexlist = indexes_choose(a);
    new_index = [1:height(data_matrix_TDS)];
    new_index(indexlist) = [];
    new_Y  = Y(new_index);
    new_X = X(new_index,:);
    indexes_choose(a)= [];
    end
    [new_XLoading,new_YLoading,new_XScore,new_YScore,new_BETA,new_PCTVAR,new_MSE,new_stats] = plsregress(new_X,new_Y,optimized_ncomp);
    new_yfit(indexlist) = [ones(size(X(indexlist,:),1),1) X(indexlist,:)]*new_BETA;
end
hold off;
figure(2); % Visualize correlation between measured and predicted responses (through leave-one-out cross-validation)
hold on;
row = 0;
for row = 1:length(Y)
    plot(Y(row),new_yfit(row),'linestyle','none','marker','o','MarkerSize',7, MarkerFaceColor='#D95319');
end
 plot([-5 5],[-5 5],'linestyle','-','linewidth',1,'color','k');
set(gca,'fontsize',20);
title('PLSR: 10-fold cross-validation');
xlabel('Measured TDS');
ylabel('Predicted TDS');
hold off;

%% Calculate VIP (Variable Importance in Projection) Scores
sum1 = zeros(1,92);
sum2 = 0;
clear SS Wnorm2;
for i = 1:optimized_ncomp
    SS(i) = (YLoading(i)^2)*(XScore(:,i)'*XScore(:,i));
end
for i = 1:optimized_ncomp
    sum2 = sum2 + SS(i);
    Wnorm2(i) = stats.W(:,i)'*stats.W(:,i);
end

clear VIP;
for counter = 1:92
    for k = 1:optimized_ncomp
        sum1(counter) = sum1(counter) + SS(k)*stats.W(counter,k)^2/Wnorm2(k);
    end
    VIP(counter) = (92*sum1(counter)/sum2)^0.5;
end

% Plot VIP scores
figure(3);
bar(VIP, 'FaceColor', '#cfa7fa');
set(gca,'XTick',1:92,'XTickLabel',vars_TDS,'YLim',[0 3.6],'Fontsize',10);
title('Variable Importance in Projection');
xlim=get(gca,'xlim');
hold on
plot(xlim,[1 1])
plot(xlim,[1.5 1.5], 'cyan')
xlabel('miRNA');
ylabel('VIP score');
xtickangle(45);
hold off;

%% table of most vip vars w scores over 1.5
vipOver1 = find(VIP(:,:)>=1.5);
varOver1 = vars_TDS(1,vipOver1);
varOver1vals = VIP(1, vipOver1);

clear varover1nicer;
varover1nicer(1, 1:3) = varOver1(1, 1:3);
varover1nicer(2, 1:3) = varOver1(1, 4:6);
varover1nicer(3, 1:3) = varOver1(1, 7:9);

figure(4)
bar(varOver1vals);
set(gca,'XTick',1:length(varOver1vals),'XTickLabel',replace(string(varOver1), "_","-"),'YLim',[0 3]);
title(['Variable with VIP scores over 1.5']);
xlabel('miRNAs');
ylabel('VIP Scores');
xtickangle(45);

%% this is the table of top 9 vip
% varover1nicer = array2table(varover1nicer);
figvip = figure;
uit = uitable(figvip,'Data',varover1nicer);
uit.ColumnName=[];
uit.RowName=[]; %removing default row numbering as in your uitable

%% write a sheet that has the top miRNAs, for later use if u want
varOver1table = cell2table(varOver1');
varOver1table.Properties.VariableNames = {'miRNAs'};
writetable(varOver1table,'miR_data/03_THCA_miR_VIPs_top_noNormal.xlsx', WriteVariableNames=true);