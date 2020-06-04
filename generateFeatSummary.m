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
featExtractTimestamp = '20200511_162714'; %'20200519_153722' (feat 3016) or '20200511_162714' (feat 3016 windows) or '20191024_122847' (feat 4548)
if strcmp(featExtractTimestamp,'20200511_162714')
    featExtractWindow = '2'; %'0','1','2'
    extractStamp = [featExtractTimestamp '_window_' featExtractWindow];
else
    extractStamp = featExtractTimestamp;
end

% load features matrix, correspondong filenames, and metadata
tierpsyFeatureTable = readtable(['/Users/sding/OneDrive - Imperial College London/aggScreening/source/features_summary_tierpsy_plate_' extractStamp '.csv'],'Delimiter',',');%,'preserveVariableNames',true);
tierpsyFileTable = readtable(['/Users/sding/OneDrive - Imperial College London/aggScreening/source/filenames_summary_tierpsy_plate_' extractStamp '.csv'],'Delimiter',',','CommentStyle','#');%,'preserveVariableNames',true);
metadataTable = readtable('/Users/sding/OneDrive - Imperial College London/aggScreening/source/metadata_aggregationScreening.csv','Delimiter',',');
% features from '20200519_153722' and ''20200511_162714' are missing the 9 food region features that need appending
if strcmp(featExtractTimestamp,'20200519_153722')
    tierpsyFeatureTable2 = readtable('/Users/sding/OneDrive - Imperial College London/aggScreening/source/features_summary_select_by_keywords_tierpsy_plate_20200526_194039_window_3.csv','Delimiter',',');%,'preserveVariableNames',true);
    tierpsyFileTable2 = readtable('/Users/sding/OneDrive - Imperial College London/aggScreening/source/filenames_summary_select_by_keywords_tierpsy_plate_20200526_194039_window_3.csv','Delimiter',',','CommentStyle','#');%,'preserveVariableNames',true);
elseif strcmp(featExtractTimestamp,'20200511_162714')
    tierpsyFeatureTable2 = readtable(['/Users/sding/OneDrive - Imperial College London/aggScreening/source/features_summary_select_by_keywords_tierpsy_plate_20200526_194039_window_' featExtractWindow '.csv'],'Delimiter',',');%,'preserveVariableNames',true);
    tierpsyFileTable2 = readtable(['/Users/sding/OneDrive - Imperial College London/aggScreening/source/filenames_summary_select_by_keywords_tierpsy_plate_20200526_194039_window_' featExtractWindow '.csv'],'Delimiter',',','CommentStyle','#');%,'preserveVariableNames',true);
end
% rename metadata column heads to match Tierpsy output
metadataTable.Properties.VariableNames{'basename'} = 'filename';

%% Join tables

% join the Tierpsy tables to match filenames with file_id. Required in case
% features were not extracted for any files.
combinedTierpsyTable = outerjoin(tierpsyFileTable,tierpsyFeatureTable,'MergeKeys',true);

% append the 9 food-region related features as necessary
if strcmp(featExtractTimestamp,'20200519_153722') | strcmp(featExtractTimestamp, '20200511_162714')
    combinedTierpsyTable2 = outerjoin(tierpsyFileTable2,tierpsyFeatureTable2,'MergeKeys',true);
    combinedTierpsyTable = outerjoin(combinedTierpsyTable,combinedTierpsyTable2,'MergeKeys',true);
end

% rename variable name to match that of metadata table if using features files from older Tierpsy versions
if strcmp(featExtractTimestamp,'20191024_122847') |  strcmp(featExtractTimestamp,'20181203_141111')
    combinedTierpsyTable.Properties.VariableNames{'file_name'} = 'filename';
end

% get just the filenames from the full path in the combined Tierpsy table
[~, fileNamesTierpsy] = cellfun(@fileparts, combinedTierpsyTable.filename, 'UniformOutput', false);
combinedTierpsyTable.filename = strrep(fileNamesTierpsy,'_featuresN','.hdf5');

% finally, join tables to get strain names for each set of features
featureTable = outerjoin(metadataTable,combinedTierpsyTable,'MergeKeys',true);

% get row logical index for valid files
rowLogInd = ~isnan(featureTable.wormNum) & featureTable.is_bad == 0 & strcmp(featureTable.is_good,'True');

% trim featureTable down to those with valid files
featureTable = featureTable(rowLogInd,:);

% export full features table
writetable(featureTable,['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fullFeaturesTable_' extractStamp '.csv']);

%% Get five worm features table with features

% trim the featureTable for 5 worm files
fiveWormLogInd = featureTable.wormNum==5;
fiveWormFeatureTable = featureTable(fiveWormLogInd,:);

% sort featureTable by strain name
strainNameColIdx = find(strcmp(fiveWormFeatureTable.Properties.VariableNames,'strain_name'));
fiveWormFeatureTable = sortrows(fiveWormFeatureTable,strainNameColIdx);

% export features table
writetable(fiveWormFeatureTable,['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fiveWorm/fiveWormFeaturesTable_' extractStamp '.csv']);

%% Get forty worm features table with features

% trim the featureTable for 40 worm files
fortyWormLogInd = featureTable.wormNum==40;
fortyWormFeatureTable = featureTable(fortyWormLogInd,:);

% sort featureTable by strain name
strainNameColIdx = find(strcmp(fortyWormFeatureTable.Properties.VariableNames,'strain_name'));
fortyWormFeatureTable = sortrows(fortyWormFeatureTable,strainNameColIdx);

% export features table
writetable(fortyWormFeatureTable,['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fortyWorm/fortyWormFeaturesTable_sw_' extractStamp '.csv']);

%% Get forty worm features table but remove all features

% remove features columns
fortyWormFeatureTable = fortyWormFeatureTable(:,1:17);

% export features table
writetable(fortyWormFeatureTable,['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fortyWorm/fortyWormFeaturesTable_' extractStamp '.csv']);