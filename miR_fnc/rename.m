%inputs -
    %wanted Names: table that has the current names in column 1 and
    %in column 2 has the names you want 
    % so example row: [ old_name   new_name]
    %data_matrix: data table with names that need to be replaced in col 1 

%outputs - 
    %a data table that is the same as data_matrix but with col 1 new names!

function renamed = rename(wantedNames, data_matrix)

names = wantedNames(:,1);
names = string(names);
replacements = wantedNames(:,2);


beforeNames = data_matrix(:,1);
lBefore = height(beforeNames);
indexChart = zeros(lBefore,1);



clear renamedIndex;
n = 1;
for i = 1:lBefore
    compare = string(beforeNames{i,1});
    x = strmatch(compare, names);
    if 0 == isempty(x)
        renamedIndex(n, 1) = i;
        x = x';
        w = length(x);
        renamedIndex(n, 2:2+w-1) = x;
        n= n+1;
    end
end

cleanedindex = renamedIndex;
h = length(renamedIndex);

if width(renamedIndex) > 2
    for i = 1:h
        m = h-i+1;
        if renamedIndex(m, 3) ~= 0
            cleanedindex(m,:) = [];
        end
    end
end

clear renamed;

rowID = cleanedindex(:,1);
row = 1;
for i = 1:height(data_matrix)
   if any(rowID == i, 1)
       renamed(row,:) = data_matrix(i,:);
       newNameIndex = cleanedindex(row,2);
       renamed{row,1} = {replacements(newNameIndex,1)};
       row = row+1;
   else
       continue;
   end
end
end