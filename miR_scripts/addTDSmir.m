%to add in the TDS scores for each tumor sample to end of data matrix
clear; clc;

%the data we want to add tds at end of:
data = readtable('miR_data/01_THCA_miR_transposed_log2_zscored_noNormal.xlsx');
obs = table2cell(data(:,1));
% obs = obs';

%% just So i could create a file with only the names and the TDS for each sample of interest
data_matrix = readtable('miR_data/01_THCA_miR_transposed_log2_zscored_noNormal.xlsx');
wantedNames = readcell('mmc5.xlsx','UseExcel',true,'Sheet','Normal_vs_Tumor DOWN');
raw = readtable('mmc3.xlsx', 'UseExcel', true, 'Sheet', 'THCA-TP (496) ');

name_and_TDS = raw(:, [1 125]);
names = table2cell(name_and_TDS(:,1));

name_and_TDS.Properties.RowNames = names;
writetable(name_and_TDS,'miR_data/02_TCHA_miR_name_and_TDS_PTC.xlsx');
%  

%% now we rename our names and convert the - to _ bc matlab silliness

namesTDSTable = readtable('02_TCHA_miR_name_and_TDS_PTC.xlsx');

h = height(namesTDSTable);
namesTDS = table2cell(namesTDSTable);

for i = 1:h
name_dash = string(namesTDS{i,1});
name_underscore = replace(name_dash, '-', '_');
namesTDS{i,1} = name_underscore;
end
%% put in the TDS vals!

w = width(data);
for i = 1:h
    sampleName = string(namesTDS{i,1});
    index = find(strcmp(obs, sampleName));
    if index == []
        continue;
    else
        data{index,w+1} = namesTDS{i,2};
    end
end
%% clear out the samples that had no TDS vals
w = width(data);
m= height(data);
for i = 1:m
    if isnan(data{m-i+1, w})
        data(m-i+1,:) =[];
    elseif data{m-i+1,w} == 0
        data(m-i+1,:) = [];
    end
end

data.Properties.VariableNames([w]) = {'TDS'};
%% save data in a sheet
writetable(data,'miR_data/02_THCA_miR_withTDS_noNormal.xlsx', 'WriteRowNames', true);