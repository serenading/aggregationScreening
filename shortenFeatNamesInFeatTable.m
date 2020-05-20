%% Function shortens the names of features that are VariableNames inside a table

% Input and output: featureTable that are n x numFeats

function featureTable = shortenFeatNamesInFeatTable(featureTable)

% separate features table into features matrix and feature names
featureMat = featureTable{:,:};
featNames = featureTable.Properties.VariableNames;
% shorten feature names
featNames = shortenFeatNames(featNames);
% put the table back together with shortened names
featureTable = array2table(featureMat,'VariableNames',featNames);

end