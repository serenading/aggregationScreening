%% This script uses three files to generate a features summary table for the 5 worm aggregation dataset,
% and a table for the 40 worm aggregation dataset in the same format except
% with no features (40 worm features to be added later on following
% calculations).

% 1. Manually curated metadata file; 2. Filenames file generated by Tierpsy; 3. Features summary file generated by Tierpsy.
% author: serenading. May 2020

clear 
close all

%% Import features data and combine with metadata

% set which feature extraction timestamp to use
featExtractTimestamp = '20191024_122847'; %'20191024_122847' or '20181203_141111'

% load features matrix, correspondong filenames, and metadata
tierpsyFeatureTable = readtable(['/Users/sding/Dropbox/aggScreening/source/features_summary_tierpsy_plate_' featExtractTimestamp '.csv'],'Delimiter',',','preserveVariableNames',true);
tierpsyFileTable = readtable(['/Users/sding/Dropbox/aggScreening/source/filenames_summary_tierpsy_plate_' featExtractTimestamp '.csv'],'Delimiter',',','preserveVariableNames',true);
metadataTable = readtable('/Users/sding/Dropbox/aggScreening/source/metadata_aggregationScreening.csv','Delimiter',',');

%% Join tables

% join the Tierpsy tables to match filenames with file_id. Required in case
% features were not extracted for any files.
combinedTierpsyTable = outerjoin(tierpsyFileTable,tierpsyFeatureTable,'MergeKeys',true);

% get just the filenames from the full path in the combined Tierpsy table
[~, fileNamesTierpsy] = cellfun(@fileparts, combinedTierpsyTable.file_name, 'UniformOutput', false);
combinedTierpsyTable.file_name = strrep(fileNamesTierpsy,'_featuresN','.hdf5');

% rename Tierpsy output and metadata column heads to match
combinedTierpsyTable.Properties.VariableNames{'file_name'} = 'filename';
metadataTable.Properties.VariableNames{'basename'} = 'filename';
   
% finally, join tables to get strain names for each set of features
featureTable = outerjoin(metadataTable,combinedTierpsyTable,'MergeKeys',true);

% export full features table
writetable(featureTable,['/Users/sding/Dropbox/aggScreening/results/fullFeaturesTable_' featExtractTimestamp '.csv']);

%% Get five worm features table with features

% trim the featureTable for 5 worm files
fiveWormLogInd = featureTable.wormNum==5 & featureTable.is_bad==0 & strcmp(featureTable.is_good,'True');
fiveWormFeatureTable = featureTable(fiveWormLogInd,:);

% sort featureTable by strain name
strainNameColIdx = find(strcmp(fiveWormFeatureTable.Properties.VariableNames,'strain_name'));
fiveWormFeatureTable = sortrows(fiveWormFeatureTable,strainNameColIdx);

% export features table
writetable(fivewormFeatureTable,['/Users/sding/Dropbox/aggScreening/results/fiveWorm/fiveWormFeaturesTable_' featExtractTimestamp '.csv']);

%% Get forty worm features table with features

% trim the featureTable for 40 worm files
fortyWormLogInd = featureTable.wormNum==40 & featureTable.is_bad==0 & strcmp(featureTable.is_good,'True');
fortyWormFeatureTable = featureTable(fortyWormLogInd,:);

% sort featureTable by strain name
strainNameColIdx = find(strcmp(fortyWormFeatureTable.Properties.VariableNames,'strain_name'));
fortyWormFeatureTable = sortrows(fortyWormFeatureTable,strainNameColIdx);

% export features table
writetable(fortyWormFeatureTable,['/Users/sding/Dropbox/aggScreening/results/fortyWorm/fortyWormFeaturesTable_sw_' featExtractTimestamp '.csv']);

%% Get forty worm features table but remove all features

% remove features columns
fortyWormFeatureTable = fortyWormFeatureTable(:,1:17);

% export features table
writetable(fortyWormFeatureTable,['/Users/sding/Dropbox/aggScreening/results/fortyWorm/fortyWormFeaturesTable_' featExtractTimestamp '.csv']);