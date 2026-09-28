function z = necessarySampleData(filename)

sampleData = readtable(filename);
clear z;


%for most data, columns one and two from tcga contain the download folder
%and file names, which is necessary for reading the data and appending it
%onto our data table
z(:,1) = sampleData(:,1); %download folder name
z(:,2) = sampleData(:,2); %file name

%for miRNA in TCGA-THCA, we want columns 3 and 4 as well, (you can change)
z(:,3) = sampleData(:,7); %sample name
z(:,4) = sampleData(:,8); %tumor classification 

end

  