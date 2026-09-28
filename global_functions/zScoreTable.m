function z = zScoreTable(file_or_table_name, numLabelCols)
%this function is to zScore our data
if istable(file_or_table_name)
    dataTable = file_or_table_name;
else 
    dataTable = readtable(file_or_table_name, 'ReadVariableNames', true);
end

clear z;

z(:,1:numLabelCols) = dataTable(:, 1:numLabelCols);

w = width(dataTable);
h = height(dataTable);
temp = zscore(dataTable{:,numLabelCols+1:w});

for row = 1:h
    for col = numLabelCols+1:w
        z{row,col} = temp(row,col-numLabelCols);
    end
end

z.Properties.VariableNames = dataTable.Properties.VariableNames;

% newFileName = strcat('z-scored-', filename);
% writetable(z, newFileName);
end