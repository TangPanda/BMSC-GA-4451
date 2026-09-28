%i think i need an ADD TDS function... its getting out of hand 

%to add in the TDS scores for each tumor sample to end of data matrix
function data = addTDS(table, sample_TDS)

%the data we want to add tds at end of:
data = table;
obs = string(data{:,1});
% obs = obs';

colTDS = array2table(zeros(length(obs),1) - 2000);
wide = width(data);
data(:, wide+1) = colTDS;


% and get a table with the sample names of interest and tds score!
namesTDSTable = sample_TDS;


h = height(namesTDSTable);
namesTDS = table2cell(namesTDSTable);

% now we rename our names and convert the - to _ bc matlab silliness
% 
% for i = 1:h
% name_dash = string(namesTDS{i,1});
% name_underscore = replace(name_dash, '-', '_');
% namesTDS{i,1} = name_underscore;
% end


%put in the TDS vals!

for i = 1:h
    sampleName = string(namesTDS{i,1});
    index = find(strcmp(obs, sampleName));
    if index == []
        continue;
    else
        data{index,wide+1} = namesTDS{i,2};
    end
end

% clear out the samples that had no TDS vals
w = width(data);
m= height(data);
for i = 1:m
    if isnan(data{m-i+1, w})
        data(m-i+1,:) =[];
    elseif data{m-i+1,w} == -2000
        data(m-i+1,:) = [];
    end
end

data.Properties.VariableNames([w]) = {'TDS'};
end