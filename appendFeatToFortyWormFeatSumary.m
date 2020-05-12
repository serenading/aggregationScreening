%% This function joins the nx2 feature table with the existing fortyWormFeatureTable in order to append new feature or update feature values
% nx2 feature table contains filenames to be used for joining in the column 1, and feature values in column 2;
% If feature name does not already exist in fortyWormFeatureTable, then
% function appends the new feature values. If feature name does exist, then
% function updates the feature values.

% author: serenading. May 2020

function [] = appendFeatToFortyWormFeatSumary(newFeatTable)

%% check that the new features table has the correct format
assert(strcmp(class(newFeatTable),'table'),'New features table must have be a "tale"')
assert(size(newFeatTable,2)==2, 'New features table must be nx2 in dimension')

%% load the most recent features table to append to
dirname = '/Users/sding/Dropbox/aggScreening/results/fortyWorm/';
d = dir([dirname 'fortyWormFeaturesTable*.csv']);
[~,idx] = max([d.datenum]); 
fortyWormFeatureTableName = d(idx).name;
fortyWormFeatureTable = readtable([dirname fortyWormFeatureTableName]);
fortyWormFeatureTableFeatNames = fortyWormFeatureTable.Properties.VariableNames;
assert(size(fortyWormFeatureTable,1)==1222,'There should be 1222 recording files for the forty-worm dataset')
assert(numel(unique(fortyWormFeatureTable.Properties.VariableNames))==size(fortyWormFeatureTable,2),...
    'Not all existing features in the fortyWormFeatureTable are unique')
n_oldFeatures = size(fortyWormFeatureTable,2)-17; % the first 17 columns are not features

%% check if the new feature already exists and update features table accordingly
newFeatName = newFeatTable.Properties.VariableNames(2);
% the new feature does not already exist, append
if nnz(strcmp(newFeatName,fortyWormFeatureTableFeatNames))==0 
    fortyWormFeatureTable = outerjoin(fortyWormFeatureTable,newFeatTable,'MergeKeys',true);
    n_newFeatures = size(fortyWormFeatureTable,2)-17;
    assert(n_newFeatures == n_oldFeatures +1, 'There should be 1 more feature than the original fortyWormFeaturesTable but this is not the case)
    disp([newFeatName ' is appended to fortyWormFeatureTable as a new feature.'])
else nnz(strcmp(newFeatName,fortyWormFeatureTableFeatNames))==1
    %%%%%% not sure what happens here if existing variable is found,
    %%%%%% whether "MergeKeys" automatically merges and keeps both
    %%%%%% versions of the values??? need to check
    fortyWormFeatureTable = outerjoin(fortyWormFeatureTable,newFeatTable,'MergeKeys',true);
    %%%%%% need to remove old feat values???
    n_newFeatures = size(fortyWormFeatureTable,2)-17;
    assert(n_newFeatures == n_oldFeatures, 'There should be the same number of features as the original fortyWormFeaturesTable but this is not the case)
    disp(['Existing values of ' newFeatName ' have been updated in the fortyWormFeatureTable.'])
else
    error('More than 1 existing features already have the same name as the new feature to be added. Cannot proceed')
end

%% save new table
timestamp=datetime('today','Format','yyyyMMdd');
writetable(fortyWormFeatureTable,['/Users/sding/Dropbox/aggScreening/results/fortyWorm/fortyWormFeaturesTable_' char(timestamp) '.csv']);