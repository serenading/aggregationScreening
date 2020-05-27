%% Function drops specified features from analysis

% Inputs: 
% featureTable
% cell array containing Tierpsy feature names as strings

% Output: 
% featureTable

function [featureTable,dropLogInd] = dropFeats(featureTable,feats2drop)

dropLogInd = false(1,size(featureTable,2));

for featCtr = 1:numel(feats2drop)
    feat = feats2drop{featCtr};
    featDropLogInd = contains(featureTable.Properties.VariableNames,feat);
    dropLogInd(featDropLogInd) = true;
    disp([num2str(nnz(featDropLogInd)) ' ' feat '-related features are dropped from analysis.'])
end

featureTable = featureTable(:,~dropLogInd);
end