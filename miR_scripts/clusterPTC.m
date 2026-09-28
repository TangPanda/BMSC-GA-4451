clear; clc;
miRNArpmData = organizeTCGA('gdc_sample_sheet_2023-01-04.xlsx'); %compile miRNA data into a table
%% save data sheet - go to miRNA_side folder to do this
writetable(miRNArpmData,'miR_data/01_THCA_miR_raw.xlsx'); %write it into a xlsx file
%% remove normal tissue solid samples
sampleIDs = miRNArpmData.Properties.VariableNames; %list of sampleIDs
sampledata = readtable('gdc_sample_sheet_2023-01-04.xlsx'); %sample data
sampleType = sampledata(:,[7,8]); %these two columns are sample ID and type

index_sampleType_normal = []; %empty array to store indexes of normal samples
for i = 1: height(sampleType)
    if strcmp(sampleType{i,2}, 'Solid Tissue Normal')
        index_sampleType_normal = [index_sampleType_normal;i];
    end
end

normalSamples = sampleType(index_sampleType_normal,:); %sample IDs of normal samples

% now get the column indexes of the normal samples in our table
colIDx_Normal = [];
for i = 1:height(normalSamples)
    normalSample_ID = normalSamples{i,1};
    for j = 1:width(sampleIDs)
        if strcmp(normalSample_ID{1}, sampleIDs{j})
            colIDx_Normal = [colIDx_Normal; j];
            break
        end
    end
end

% and rewrite the datatable without normal samples 
allIDx = [1:width(miRNArpmData)];
tumorIDx = setdiff(allIDx,index_sampleType_normal);

noNormalmiRNArpmData = miRNArpmData(:,tumorIDx);

%save it as a file
writetable(noNormalmiRNArpmData,'miR_data/01_THCA_miR_raw_noNormal.xlsx'); %write it into a xlsx file

%% reduce data based on TCGA-THCA paper how they renamed miRs
%now we get only the data that they used in their experiment and rename as
%they did 
[status,sheets] = xlsfinfo('mmc5.xlsx'); %from TCGA-THCA supplemental info
sheets = sheets(:,3:18);

%this is to grab the miRs TCGA-THCA classified and their 'official name'
clear wantedRename;
n = 1;
for i = 1:16
    sheetName = string(sheets{1,i});
    wantedRename1 = readcell('mmc5.xlsx','UseExcel',true,'Sheet',sheetName);
    wantedRename1 = wantedRename1(2:end, 1:2);
    h = length(wantedRename1);
    wantedRename(n:n+h-1,1:2) = wantedRename1(:,1:2);
    n = n+h;
end
wantedRename = cell2table(wantedRename);
uN = unique(wantedRename);
uN = table2cell(uN);
%% Rename our data table 
data_matrix = readtable('miR_data/01_THCA_miR_raw_noNormal.xlsx');
data_matrix = rename(uN, data_matrix);
writetable(data_matrix,'miR_data/01_THCA_miR_miRsRenamed_noNormal.xlsx');
%% transpose the data to get our correct vars and obs
%we have to transpose the data to get our correct vars and obs :-)
data_matrix = readtable('01_THCA_miR_miRsRenamed_noNormal.xlsx');
data_matrix_t = transposeTable(data_matrix);
writetable(data_matrix_t,'miR_data/01_THCA_transpose_miR_miRsRenamed_noNormal.xlsx', 'WriteRowNames', true, 'WriteVariableNames', true)


%% log2normalize and zscore

data_matrix_t = readtable('01_THCA_transpose_miR_miRsRenamed_noNormal.xlsx');

%because we have rownames, we use log2tableROWNAME function
miRNAlog2 = log2tableROWNAME(data_matrix_t);

writetable(miRNAlog2,'miR_data/01_THCA_miR_transposed_log2_noNormal.xlsx');

%% zscore the log normalized data and create a new excel sheet for it 
%zscore the log normalized data and create a new excel sheet for it :-)
miRNAlog2Zscore = zScoreTable('01_THCA_miR_transposed_log2_noNormal.xlsx', 1);
writetable(miRNAlog2Zscore,'miR_data/01_THCA_miR_transposed_log2_zscored_noNormal.xlsx');

%% establish variables, observations, and the data matrix
data_matrix_tocluster = readtable('01_THCA_miR_transposed_log2_zscored_noNormal.xlsx');
observations = data_matrix_tocluster{:,1};
vars = data_matrix_tocluster.Properties.VariableNames([2:end]);
vars = replace(string(vars), "_", "-");
% vars = vars(1,2:end);
data_matrix_vals = data_matrix_tocluster{:,2:end};
%% this part is to actually cluster our data!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cluster1 = clustergram(data_matrix_vals,'Cluster',3,'Symmetric',true,'Colormap',redblue(150),...
    'DisplayRange',3,'RowPDist','correlation','ColumnPDist','euclidean','Linkage','complete',...
    'RowLabels',observations,'ColumnLabels',vars, 'ColumnLabelsRotate', 45); 

plot(cluster1)
hold on
c = colorbar('Ticks',[-1,0,1], ...
         'TickLabels',["", "low expression","high expression"]);
% c.Limits = [-1,1];
%% random stuff with clusters and dendrograms (not relevant)
% set(cluster1,'Dendrogram',1.5);
% rm = struct('GroupNumber',{383,384,385,380, 368},'Annotation',{'Group383','Group384','Grouop385','Group380','Group368'},...
%      'Color',{'g','y', 'b','m','r'});
% set(cluster1, 'RowGroupMarker', rm);
% %% i dunno how to automate this but - get the groups 560, 566, 567, 569
% 
% group383 = get(Group383, "RowLabels");
% group384 = get(Group384, "RowLabels");
% group385 = get(Group385, "RowLabels");
% % group386 = get(Group386, "RowLabels");
% group380 = get(Group380, "RowLabels");
% group368 = get(Group368, "RowLabels");
% 
% %%
% fiveGroups = string(group383);
% fiveGroups(1:length(group384),2) = string(group384);
% fiveGroups(1:length(group385), 3) = string(group385);
% fiveGroups(1:length(group380),4) = string(group380);
% fiveGroups(1:length(group368),5) = string(group368);
% 
% fiveGroups = array2table(fiveGroups);
% %%
% row_labels_in_order = cluster1.RowLabels;
% TDS_no_order = data_matrix(:, [1 94]);
% 
% %rearrange the TDS scores
% TDS_rearranged = rearrange(TDS_no_order,row_labels_in_order);
% TDS_rearranged = flip(TDS_rearranged);
% TDS_rearranged_vals = table2array(TDS_rearranged(:,2));
% 
% heatmaparrayTDS = TDS_rearranged_vals;
% xvalues = 'TDS';
% yvalues = string(TDS_rearranged{:,1});
% 
% figure(2)
% htds = heatmap(xvalues,yvalues,heatmaparrayTDS, 'CellLabelColor','none');
% 
% htds.Title = 'TDS';
% htds.XLabel = 'score';
% htds.YLabel = 'Cell ID';
% htds.Colormap = redblue;
% clim([-4.5 4.5])
% 
% %%
% fiveGroups.Properties.VariableNames = [{'Group383'} {'Group384'} {'Group385'} {'Group380'} {'Group368'}];
% writetable(fiveGroups, 'fiveGroups_clustergram.xlsx', 'WriteVariableNames',true);


%%
    function insertColorbarCBALWAYS(obj)

         hFig= gcbf;

         obj.Colorbar = true;

    end