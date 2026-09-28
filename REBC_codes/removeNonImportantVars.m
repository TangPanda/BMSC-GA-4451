%function to remove nonimportant things
function important = removeNonImportantVars(tableWithVars, notImportantVars)

important = tableWithVars;
clear notimpIDx;
notimpIDx = [];
for i = 1:length(notImportantVars)
    not_imp_name = notImportantVars(i);
    idx = find(strcmp(string(important.Properties.VariableNames'),not_imp_name));
    if ~isnan(idx)
        notimpIDx = [notimpIDx; idx];
    end
end

notimpIDx = sort(notimpIDx);

%% now remove those indices
w = width(important);
for m = 1:w    
    if ismember(w-m+1, notimpIDx)
        important(:, w-m+1) = [];
    end
end
end
