%function to transpose data tables - not very useful because its very
%straightforward, but it cluttered up the code in the main script haha

function transposed_data_matrix = transposeTable(table)
data_matrix = table;

%observations are in column 1
obs = table2array(data_matrix(:,1));

%convert data to array
arrayData = table2array(data_matrix(:,2:end));
%transpose
transposed_data_matrix = array2table(arrayData');

%row names are now var names except not from col 1
transposed_data_matrix.Properties.RowNames = data_matrix.Properties.VariableNames(2:end);
%var names are thhe observations from column 1 of old table
obs_transpose = obs';
transposed_data_matrix.Properties.VariableNames = string(obs_transpose);

%return transposed data table with correct vars and obs labels
end
