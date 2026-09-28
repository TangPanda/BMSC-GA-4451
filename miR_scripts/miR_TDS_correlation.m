clear; clc;


%Thinking out loud:
% want to correlate TDS with miRNAs that are most important... how to get
%most important miRNAs?
%maybe penalized regressive models? or we can start off with just the top 9
%miRNAs that we identified in PLSR! Let's start there

%spearman correlate for miRNAS we got from PLSR
miRNAs = readtable('miR_data/miR_VIPs_top.xlsx', 'ReadVariableNames',true);
h = height(miRNAs);

%general table with the sample name and the TDS in cols 1 and 2, resp:
allmiR_data_log2zscoreTDS = readtable('miR_data/TCGAmiRNAwithTDS_PTC.xlsx');
obsTDS_table = cell2table(allmiR_data_log2zscoreTDS.obs);
obsTDS_table(:,2) = array2table(allmiR_data_log2zscoreTDS.TDS);
obsTDS_table.Properties.VariableNames = [{'obs'} {'TDS'}];

%%
%make a correlation for each miRNA and a table for each!
clear corrtable;
corrtable = array2table(obsTDS_table.TDS);
corrtable.Properties.VariableNames = {'TDS'};
corrtable.Properties.RowNames = obsTDS_table{:,1};
for i = 1:h
    miRNA = string(miRNAs{i,1});
    corrtable{:,i+1} = allmiR_data_log2zscoreTDS{:, miRNA};
    varname_miRNA = replace(miRNA, "miR_", "");
    varname_miRNA = char(replace(varname_miRNA,"_","-"));
    corrtable.Properties.VariableNames(i+1) = {varname_miRNA};
end


%%
figure(4)
[r, pvalue, handle] = corrplot(corrtable, 'varName', corrtable.Properties.VariableNames, ...
    Type = 'Spearman', TestR = "on");
p = pvalue(1,:);
%% 

fig = uifigure;
figure(5)
uit = uitable(fig,'Data',p);
% loop to only show first row of corrplot but im not sure how to fully
% clean
% for i = 2:10
%     for j = 1:10
%         delete(subplot(10, 10, (i-1)*10+j));
%         gca()
%     end
% end
%%  fancy guy
%i wanna make one of those cool nice fancy circle dot fancy guys

TDS_corrdata = r{1,:};
clrLim = [-1,1];
diamLim = [.1, 1];
labels = corrtable.Properties.VariableNames;

x = 1 : 1 : size(TDS_corrdata, 2);
[xAll, yAll] = meshgrid(x,1);
xAll(TDS_corrdata == 0) = nan;

cmap = redblue(256);
corrdata_scaled = (TDS_corrdata - clrLim(1))/range(clrLim);
clrID = discretize(corrdata_scaled, linspace(0,1,size(cmap,1)));

corrdata_scaled = (abs(TDS_corrdata));
diamSize = corrdata_scaled * range(diamLim) + diamLim(1);
diamSize(1,1) = 0;

fancy_fig_1 = figure();
ax = axes(fancy_fig_1);
hold(ax, 'on')
colormap(redblue);
tickvalues = 1:length(TDS_corrdata);
% z = zeros(size(tickvalues));
text(tickvalues, yAll+.5, labels, 'HorizontalAlignment', 'right','Rotation',45);

theta = linspace(0,2*pi,50); % the smaller, the less memory req'd.
h = arrayfun(@(i)fill(diamSize(i)/2 * cos(theta) + xAll(i), ...
    diamSize(i)/2 * sin(theta) + yAll(i), cmap(clrID(i),:),'LineStyle','none'),1:numel(xAll));
axis(ax,'equal')
axis(ax,'tight')
set(ax,'YDir','Reverse')
title('Correlation between TDS and miRNAs')
colorbar()
caxis(clrLim);
axis off
%% what if i did corr fancy guy for all the mirnas vs ea other :-O?
all_corrdata = tril(table2array(r),-1);
all_corrdata(logical(eye(size(all_corrdata)))) = 1;

x = 1 : 1 : size(all_corrdata, 2);
y = 1 : 1 : size(all_corrdata,1);
[xAll, yAll] = meshgrid(x,y);
xAll(all_corrdata == 0) = nan;

cmap = redblue(256);
all_corrdata_scaled = (all_corrdata - clrLim(1))/range(clrLim);
clrID = discretize(all_corrdata_scaled, linspace(0,1,size(cmap,1)));

all_corrdata_scaled = (abs(all_corrdata));
diamSize = all_corrdata_scaled * range(diamLim) + diamLim(1);

fancy_fig = figure();
ax = axes(fancy_fig);
hold(ax, 'on')
colormap(redblue);
tickvalues = 1:length(all_corrdata);
% z = zeros(size(tickvalues));
xtick = zeros(1, length(x));
ytick = zeros(1, length(y)) + length(tickvalues) + .5;
text(tickvalues, ytick, labels, 'HorizontalAlignment', 'right','Rotation',45);
text(xtick, tickvalues, labels, 'HorizontalAlignment', 'right','Rotation',45);

theta = linspace(0,2*pi,50); % the smaller, the less memory req'd.
h = arrayfun(@(i)fill(diamSize(i)/2 * cos(theta) + xAll(i), ...
    diamSize(i)/2 * sin(theta) + yAll(i), cmap(clrID(i),:),'LineStyle','none'),1:numel(xAll));
axis(ax,'equal')
axis(ax,'tight')
set(ax,'YDir','Reverse')
title('Pairwise Correlation between miRNAs and TDS')
colorbar()
caxis(clrLim);
axis off


