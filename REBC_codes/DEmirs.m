% differential expression 
clear; clc;

tumorXP= readtable("REBC_test/REBC_data/mirtumor_log2.csv", 'ReadVariableNames',true);
normalXP = readtable("REBC_test/REBC_data/mirnorm_log2.csv", 'ReadVariableNames',true);
%%
t_tumorXP = transposeTable(tumorXP);
t_normalXP = transposeTable(normalXP);
t_tumorXP_data = table2array(t_tumorXP);
t_normalXP_data = table2array(t_normalXP);
%% differential expression 

pvals = mattest(t_tumorXP_data, t_normalXP_data);
t_tumorDE = t_tumorXP(find(pvals<=0.05),:);
t_normDE = t_normalXP(find(pvals<=0.05),:);
% %% write sheet
% writetable(t_tumorDE, "REBC_test/REBC_data/tDEmirs_tumor.csv", "WriteVariableNames",true, "WriteRowNames",true);
% writetable(t_normDE, "REBC_test/REBC_data/tDEmirs_normal.csv", "WriteVariableNames",true,"WriteRowNames",true);

%% data org
clear; clc;
t_tumorDE = readtable("REBC_test/REBC_data/tDEmirs_tumor.csv", "ReadVariableNames",true);
t_normDE = readtable("REBC_test/REBC_data/tDEmirs_normal.csv", "ReadVariableNames",true);

tumorDE = transposeTable(t_tumorDE);
normalDE = transposeTable(t_normDE);
%% zscore and write into sheets
zL2_tumormirs = zScoreTable(tumorDE, 1);
zL2_normmirs = zScoreTable(normalDE, 1);

writetable(zL2_tumormirs, "REBC_test/REBC_data/DEmirtumor_l2zscored.csv", "WriteVariableNames",true, "WriteRowNames", true);
writetable(zL2_normmirs, "REBC_test/REBC_data/DEmirnorm_l2zscored.csv", "WriteVariableNames",true,  "WriteRowNames", true);
%%  rename and add TDS 
clear zL2_normmirs zL2_tumormirs;
zL2_tumormirs = readtable("DEmirtumor_l2zscored.csv", "ReadVariableNames", true);
zL2_normmirs = readtable("DEmirnorm_l2zscored.csv", "ReadVariableNames", true);
mirTableData = readtable('REBCsample_and_mir_tabledata.xlsx', "ReadVariableNames", true);

wantedNames = string([mirTableData.Var5, mirTableData.Var4]);

renamed_tumormirs = rename(wantedNames, zL2_tumormirs);
renamed_normmirs = rename(wantedNames, zL2_normmirs);

% %% write renamed into sheets
writetable(renamed_tumormirs, "REBC_test/REBC_data/DErenamed_mirtumor_l2zscored.csv", "WriteVariableNames",true);
writetable(renamed_normmirs, "REBC_test/REBC_data/DErenamed_mirnorm_l2zscored.csv", "WriteVariableNames",true);

%%
namesTDSTable = readtable("REBC_test/REBC_data/REBCsample_and_TDS.xlsx");
tumorwithTDS = addTDS(renamed_tumormirs, namesTDSTable);
normwithTDS = addTDS(renamed_normmirs, namesTDSTable);

% %% write into sheets
writetable(tumorwithTDS,"REBC_test/REBC_data/DEmirREBC_tumor_TDS.csv", WriteVariableNames=true);

writetable(normwithTDS,"REBC_test/REBC_data/DEmirREBC_norm_TDS.csv", WriteVariableNames=true);


