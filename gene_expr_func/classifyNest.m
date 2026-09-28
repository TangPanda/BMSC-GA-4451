%%
function classificationTable = classifyNest(tbl, vars, rows)
clear classificationTable;
for i = 1:length(rows)
    row = rows(i);
    subrows = unique(tbl.(row));
    idx = find(cellfun(@(s) any(ismember(['X', '-'], s)), subrows));
    subrows(idx) = [];
    for j = 1:length(vars)
        var = vars{j};
        lowJ = tbl.(var) < 0;
        highJ = tbl.(var) >= 0; 
        clear subTable
        for k = 1:length(subrows)
            clear subRowYN;
            subrow = subrows{k};
%             if ismember('X', subrow)
%                 continue;
%             elseif ismember('-', subrow)
%                 continue;
%             else
%                 subRowYN = matches(tbl.(row), subrow);
%                 subTable(k,1) = sum(lowJ & subRowYN);
%                 subTable(k,2) = sum(highJ & subRowYN);
%             end
             subRowYN = matches(tbl.(row), subrow);
             subTable(k,1) = sum(lowJ & subRowYN);
             subTable(k,2) = sum(highJ & subRowYN);
        end
        subTable = array2table(subTable);
        subTable.Properties.VariableNames = [{char("l_" + var)} {char("h_" +var)}];
        subTable.Properties.RowNames = subrows;
        classificationTable(i,j) = {subTable};
    end
end
classificationTable = cell2table(classificationTable);
classificationTable.Properties.RowNames = rows;
classificationTable.Properties.VariableNames = vars;
end