% this function will help you organize TCGA data

%inputs - sample sheet - file roster that contains download data info
%outputs - compiled table of data from TCGA with row and var names
function organizedTable = organizeTCGA(sampleDataFileName)

%example datasheet: necessarySampleData('gdc_sample_sheet.xlsx');


%read the datasheet that has your sample info/data in it

%go to necessarySampleData function to change which columns u want to keep
%from sample sheet (what info is relevant)
datasheet = necessarySampleData(sampleDataFileName);

%choose a random file to get the column of miRNAs and also dimentions:
fname1 = string(datasheet{1,2});
s = readtable(fname1);
h = height(s);
w = height(datasheet) - 1;

%put that column of miRNAs into our organizedTable
clear organizedTable;
organizedTable(:,1) = s(:,1);

%loop through all our data files to get the rpm of each miRNA
for id = 1:w
    fname = string(datasheet{id,2});
    info = readtable(fname);
    sampleID = string(datasheet{id,3});
    organizedTable(:, id+1) = info(:,3);
    organizedTable.Properties.VariableNames(id+1) = [sampleID];

end

end
    
    
    
