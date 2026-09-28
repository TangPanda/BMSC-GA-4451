function bestalphaTable  = findBestAlpha(x,y,vars,nLoops,stepsize,kfold)
clear MSEtable;
MSEtable = zeros(1/stepsize,2);
stepsarray = stepsize:stepsize:1;

n = length(stepsarray);
parfor i = 1:n
    alpha = stepsarray(i);
    bcoeftable50 = loopRegressCoefs(x,y,vars,nLoops,alpha,kfold);
    MSEtable(i,:) = [alpha, mean(bcoeftable50{height(bcoeftable50),:})];
    disp([string(i) "out of " string(n)]);
end
%%
% bestalphaTable = MSEtable(find(min(MSEtable(:,2))),1);
bestalphaTable = MSEtable;

end