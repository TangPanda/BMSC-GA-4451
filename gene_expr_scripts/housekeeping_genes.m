 clear; clc;

% these are the housekeeping genes we know of: 
% beta-actin (actb), gapdh, hprt1, tbp, gusb, HMBS, TFRC,SDHA,YWHAZ,B2M,
% ARBP,RPLP0,RPL13A,PPIA,PGK1,RRN18S



sharedgenes = readtable('shared_genes.csv', 'ReadVariableNames',false);
sharedgenes = string(table2array(sharedgenes));
housekeeping = ["ACTB"; "GAPDH"; "HPRT1"; "TBP"; "GUSB"; "HMBS"; "TFRC"; ...
    "B2M"; "ARBP"; "RPLP0"; "RPL13A"; "PPIA"; "PGK1"; "RRN13S"; "SDHA"; "YWHAZ"];
ind = zeros(length(housekeeping),1);
for i = 1:length(housekeeping)
    hkg = housekeeping(i,1);
    ID = find(sharedgenes == hkg);
    if isempty(ID)
        continue;
    else
        ind(i,1) = ID;
    end
end

