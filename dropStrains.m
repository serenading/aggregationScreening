%% Function drops specified strains from analysis

% Inputs: 
% featureTable
% cell array containing strain names as strings

% Output: 
% featureTable

function [featureTable,dropLogInd] = dropStrains(featureTable,strains2drop)

dropLogInd = false(1,size(featureTable,1));

for strainCtr = 1:numel(strains2drop)
    strain = strains2drop{strainCtr};
    strainDropLogInd = strcmp(featureTable.strain_name,strain);
    dropLogInd(strainDropLogInd) = true;
    disp([strain ' is dropped from analysis.'])
end

featureTable = featureTable(~dropLogInd,:);
end