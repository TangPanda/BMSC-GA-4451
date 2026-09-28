function t = log2table(dataTable)

t = dataTable(:,1);
dataTablevals = table2array(dataTable(:,2:end));
h = height(dataTablevals);
w = width(dataTablevals);
for i = 1:h
    for j = 1:w
        data = dataTablevals(i,j);
        t{i,j+1} = log2(data+1);
    end
end
% t = cell2table(t);
vars = dataTable.Properties.VariableNames;
% vars{1,1} = 'obs';
t.Properties.VariableNames = vars;
t.Properties.RowNames = dataTable.Properties.RowNames;

end