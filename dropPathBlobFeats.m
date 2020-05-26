%% Function drops path-related and blob features from analysis

% Input & Output: cell array containing Tierpsy feature names

function featureTable = dropPathBlobFeats(featureTable)

featLogInd = ~contains(featureTable.Properties.VariableNames,'path') &...
    ~contains(featureTable.Properties.VariableNames,'blob');
featureTable = featureTable(:,featLogInd);

end