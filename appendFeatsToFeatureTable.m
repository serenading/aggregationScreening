%% This function appends newFeatureTable with the existing featureTable.
% This function currently only works with 40 worms, but can easily be modified to work with 5 worms. 
% Duplicate features are removed from newFeatureTable before appending. 
% Tables are joined using filename field.

% author: serenading. May 2020

function [] = appendFeatsToFeatureTable(newFeatureTable,wormNum,extractStamp)

%% Check that the new features table has the correct format
assert(isa(newFeatureTable,'table'),'New feature table must have be in table format.')
assert(ismember('filename',newFeatureTable.Properties.VariableNames), 'New feature table must contain filenames column to enable table joining.')

%% Load the most recent features table to append to
if wormNum == 40
    dirname = '/Users/sding/OneDrive - Imperial College London/aggScreening/results/fortyWorm/';
    d = dir([dirname 'fortyWormFeaturesTable_' extractStamp '*.csv']); 
elseif wormNum ==5
    dirname = '/Users/sding/OneDrive - Imperial College London/aggScreening/results/fiveWorm/';
    d = dir([dirname 'fiveWormFeaturesTable_' extractStamp '*.csv']); 
end

[~,idx] = max([d.datenum]); % this gets the index of the latest file
featureTableName = d(idx).name;
featureTable = readtable([dirname featureTableName]);
n_oldFeatures = size(featureTable,2)-17;

%% Check that the existing table is in the correct format
featNames = featureTable.Properties.VariableNames;
if wormNum ==40
    assert(size(featureTable,1)==1222,'There should be 1222 recording files for the forty-worm dataset')
    assert(numel(unique(featureTable.filename))==1222,'There should be 1222 recording files for the forty-worm dataset')
elseif wormNum ==5
    assert(size(featureTable,1)==787,'There should be 787 recording files for the forty-worm dataset')
    assert(numel(unique(featureTable.filename))==787,'There should be 787 recording files for the forty-worm dataset')
end
assert(numel(unique(featNames))==size(featureTable,2),...
    'Some existing features in the old featureTable are duplicated')

%% Remove duplicate features from the new features table (necessary for outerjoin without creating extra feature columns)
newFeatureTableFilenames = newFeatureTable.filename; % save filename to append for joining
newFeatNames = newFeatureTable.Properties.VariableNames;
both = ismember(newFeatNames,featNames); % get indices for duplicate features
newFeatureTable = newFeatureTable(:,~both); % remove duplicate features
newFeatureTable.filename = newFeatureTableFilenames; % re-append filenames for joining

%% Join old and new feature tables
featureTable = outerjoin(featureTable,newFeatureTable,'MergeKeys',true);
n_newFeatures = size(featureTable,2)-17;
n_featDiff = n_newFeatures - n_oldFeatures;
assert(n_featDiff == size(newFeatureTable,2)-1,'Incorrect number of features appended.')
disp([num2str(n_featDiff) ' new features are appended to featureTable. n_oldFeatures: ' num2str(n_oldFeatures) '; n_newFeatures: ' num2str(n_newFeatures) '.'])

%% save new table as a .csv file
timestamp=datetime('today','Format','yyyyMMdd');
if wormNum == 40
    savename = [dirname 'fortyWormFeaturesTable_' extractStamp '_' char(timestamp) '_ft' num2str(n_newFeatures) '.csv'];
elseif wormNum == 5
    savename = [dirname 'fiveWormFeaturesTable_' extractStamp '_' char(timestamp) '_ft' num2str(n_newFeatures) '.csv'];
end
writetable(featureTable,savename)