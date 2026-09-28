%this function is for gene expression data tables that have row names

%input - dataTable with row names

%output - dataTable that has been log2'd , with row names in column 1

function t = log2tableROWNAME(dataTable)
h = height(dataTable);
w = width(dataTable);
clear t;
maybeRowNames = dataTable.Properties.RowNames;

%this if statement is kinda useless bc we know maybeRowNames is not empty
if isempty(maybeRowNames) == 1
    t(:,1) = table2array(dataTable(:,1));
else
    t(:,1) = maybeRowNames(:,1);
end

%log2(data+1) - the standard for gene expression (TPM)
for i = 1:h
    for j = 2:w
        data = dataTable{i,j};
        t{i,j} = log2(data+1);
    end
end

%write new table
t = cell2table(t);
vars = dataTable.Properties.VariableNames;
vars{1,1} = 'obs';
t.Properties.VariableNames = vars;

end