%% loop it 10, 50, or 100 times!
function coefficient_table = loopRegressCoefs(x,y,vars,nLoops, alpha, k_cv)

clear coefficient_table;

for i = 1:nLoops

    [B,FitInfo] = lasso(x,y,'CV',k_cv, 'PredictorNames', vars, 'Alpha', alpha);
    
    blam = FitInfo.Index1SE;
    bcoefs = B(:,blam);
    coefficient_table(1:length(bcoefs)+1,i) = [bcoefs; FitInfo.MSE(blam)];

end

coefficient_table = array2table(coefficient_table);
coefficient_table.Properties.RowNames = [vars {'MSE'}];
end





