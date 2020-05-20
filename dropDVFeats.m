%% Function drops dorsal ventral features from analysis

% Input & Output: cell array containing Tierpsy feature names

function featuresTable = dropDVFeats(featuresTable)

% get a list of valid feature names
features = readtable('auxiliary/ft_names_without_ventrally_signed.csv');
featureNames = features.Variables;
% shorten feature names
featureNames = shortenFeatNames(featureNames);
% trim down features matrix
featureMatrix = featureTable{:,featureNames};
% put the matrix back into a table
featureTable = array2table(featureMatrix,'VariableNames',featureNames);

end