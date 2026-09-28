function z = necessarySampleData(filename)

sampleData = readtable(filename);
clear z;


%for most data, columns one and two from tcga contain the download folder
%and file names, which is necessary for reading the data and appending it
%onto our data table
z(:,1) = sampleData(:,1); %download folder name
z(:,2) = sampleData(:,2); %file name

%for rna seq in TCGA-THCA, we want columns 7 and 8 as well, (you can change)
z(:,3) = sampleData(:,4); %data type
z(:,4) = sampleData(:,6); %Case ID
z(:,5) = sampleData(:,7); %sample name
z(:,6) = sampleData(:,8); %tumor classification 

end

  