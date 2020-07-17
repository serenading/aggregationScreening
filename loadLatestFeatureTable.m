%% Functin loads the latest featureTable, because these get updated all the time as new features are added.
% Author: @serenading. July 2020

function [featureTable, pathname] = loadLatestFeatureTable(extractStamp,wormNum)

if wormNum == 40
    dirname = '/Users/sding/OneDrive - Imperial College London/aggScreening/results/fortyWorm/';
    d = dir([dirname 'fortyWormFeaturesTable_' extractStamp '*.csv']); 
elseif wormNum ==5
    dirname = '/Users/sding/OneDrive - Imperial College London/aggScreening/results/fiveWorm/';
    d = dir([dirname 'fiveWormFeaturesTable_' extractStamp '*.csv']);
else
    error('Please specify a valid wormNum.')
end

[~,idx] = max([d.datenum]); % this gets the index of the latest file
featureTableName = d(idx).name;
pathname = [dirname featureTableName];
featureTable = readtable(pathname,'Delimiter',',','preserveVariableNames',true);

end