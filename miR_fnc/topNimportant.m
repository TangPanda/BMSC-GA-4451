%input - 
    %table: a data table you want to get the top important vars of
    %uptocolumn: we will get top imp vars of columns 1 to (uptocolumn)
    %n: the number of most important vars you want

%output - topNtable: table of the indices where the top N lie and their vals, in order!

function topNtable = topNimportant(table, uptocolumn, n)

    clear topNtable
    for i = 1:uptocolumn
        loadings = abs(table(:,i));
        indexColumn = (i-1)*2+1;
        [loadVal,top5] = maxk(loadings, n);
        topNtable(:, indexColumn) = top5;
        loadingColumn = i*2;
        topNtable(:, loadingColumn) = loadVal;
    end
    topNtable = array2table(topNtable);

end