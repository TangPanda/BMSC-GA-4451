function markerColors = addTDSColor(data_matrix, n)
cmap = colormap(redblue(n));
yourValues = data_matrix;
v = rescale(yourValues, 1, n); % Nifty trick!
numValues = length(yourValues);
markerColors = zeros(numValues, 3);
% Now assign marker colors according to the value of the data.
for k = 1 : numValues
    row = round(v(k));
    markerColors(k, :) = cmap(row, :);
end


end