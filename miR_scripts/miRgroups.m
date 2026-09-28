clear; clc;

%we r going to find the miRNA data for the four groups found in clustergram

%read in the groups!
groups = readtable('fourGroups_clustergram.xlsx', ReadVariableNames=true);

miRNA_data_all = readtable('TCGAmiRNAwithTDS_PTC.xlsx');

%% make tables for the groups with their miRNAs and TDS
%and now make tables for the groups:
% group560miR = getSpecificVals(groups(:,1), miRNA_data_all);

group383 = string(groups{:,1});
group383_dub = [group383 group383];

group383_miR = rename(group383_dub, miRNA_data_all);

group384_miR = rename([string(groups{:,2}) string(groups{:,2})], miRNA_data_all);
group385_miR = rename([string(groups{:,3}) string(groups{:,3})], miRNA_data_all);
% group380_miR = rename([string(groups{:,4}) string(groups{:,4})], miRNA_data_all);
% group368_miR = rename([string(groups{:,5}) string(groups{:,5})], miRNA_data_all);
group386_miR = rename([string(groups{:,4}) string(groups{:,4})], miRNA_data_all);


%% and heatmap thhe TDS for the groups
heatmaparray = table2array(group386_miR(:,94));
xvalues = group386_miR.Properties.VariableNames(94);
yvalues = string(group386_miR{:,1});

figure(1)
hc = heatmap(xvalues,yvalues,heatmaparray);

hc.Title = 'TDS group 386';
hc.XLabel = 'TDS';
hc.YLabel = '';
hc.Colormap = redblue;
clim([-4.5 4.5])
Ax = gca;
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));

%% swarm chart for the top 9 miRNAs from PLSR

miRNAs = string(table2array(readtable('miR_VIPs_top.xlsx')));

%% box and whisker plot

figure(2)
h = length(miRNAs);
for i = 1:h
    subplot(3,3,i);
    hold on;
    clear miRNA_group;
    miR = miRNAs(i,1);
    miRNA_group = group383_miR(:,miR);
    miRNA_group(1:height(group384_miR),2) = group384_miR(:,miR);
    miRNA_group(1:height(group385_miR),3) = group385_miR(:,miR);
    miRNA_group(1:height(group386_miR),4) = group386_miR(:,miR);
    % miRNA_group(1:height(group368_miR),5) = group368_miR(:,miR);
    
    miRNA_group = table2array(miRNA_group);
    miRNA_group(miRNA_group == 0) = nan;
    
    boxplot(miRNA_group, ['Group383'; 'Group384'; 'Group385'; 'Group386'], ...
        'Notch','on')
    title(['miRNA expression distribution of ' replace(miR,"_","-")]);
    xlabel('group');
    ylabel(replace(miR,"_","-"));
    ylim([-3.25 3.25]);
    legend('off');
end
%% kruskal-wallis test 

h = length(miRNAs);
clear pvaltable;
for i = 1:h
    clear miRNA_group;
    miR = miRNAs(i,1);
    miRNA_group = group383_miR(:,miR);
    miRNA_group(1:height(group384_miR),2) = group384_miR(:,miR);
    miRNA_group(1:height(group385_miR),3) = group385_miR(:,miR);
    miRNA_group(1:height(group386_miR),4) = group386_miR(:,miR);
    % miRNA_group(1:height(group368_miR),5) = group368_miR(:,miR);
    
    miRNA_group = table2array(miRNA_group);
    miRNA_group(miRNA_group == 0) = nan;
%     [p,tbl,stats] = anova1(miRNA_group);
    figure(i)
    miR = replace(miR, "_", '-');
    [p,tbl,stats] = kruskalwallis(miRNA_group, ["group383" "group384" "group385" "group386"],"off");
    pvaltable((i-1)*7+1,1) = array2table(miR);
    results = multcompare(stats);
    results(results==1) = 383;
    results(results==2) = 384;
    results(results==3) = 385;
    results(results==4) = 386;
    results = array2table(results(:,:));
    pvaltable((i-1)*7+2:i*7, 1:6) = table2cell(results(:,:));
    title(['comparison of ' miR ' across groups']);
%     figure(2*i)
%     uitable('Data', results,'ColumnName', {'GroupA', 'GroupB', 'Lower Limit', 'A-B', 'Upper Limit', 'p-Value'});
%     tbl = array2table(results,"VariableNames", ...
%     ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])
end

%%
pvaltable.Properties.VariableNames = [{'GroupA'} {'GroupB'} {'Lower Limit'} {'A-B'} {'Upper Limit'} {'p-Value'}];
writetable(pvaltable, 'kruskal-wallis test data.xlsx', 'WriteVariableNames',true);
