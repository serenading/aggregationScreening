clear
close all

addpath('auxiliary/')

%% Script extracts Tierpsy blob features from 40 worm tracking data to generate a feature table for joining onto the master feature table

% author: serenading. June 2020.

%% Specify analysis parameters

% which HD is plugged in right now?
HD1 = true; % true or false

% which base features to expand and add to features table?
featureBaseNames = {'blobSpeed','d_blobSpeed', 'd_blobSpeed_abs','blobArea','blobCompactness','blobPerimeter','blobSolidity','blobQuirkiness','blobHu0','blobHu1','blobHu2','blobHu3','blobHu4','blobHu5','blobHu6','blobLength','blobWidth','blobOrientation'};
n_expandedFeats = 100; % number of expanded features per base feature

% set which feature extraction timestamp to use (doesn't matter as we are calculating our own feats)
extractStamp = '20200519_153722';

% how many times the area of a single worm counts as a cluster?
clusterArea = 3;

% some file attributes
frameRate = 25; % 25 frames per second, read by frameRate = double(h5readatt(strrep(filename,'featuresN','skeletons'),'/plate_worms','expected_fps'));
pixelToMicron = 10; % 10 microns per pixel, read by pixelsize = double(h5readatt(filename,'/trajectories_data','microns_per_pixel'))

%% Load feature table and find which 40 worm files are on which HD
% load forty worm feature table
featureTable = readtable(['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fortyWorm/fortyWormFeaturesTable_' extractStamp '.csv'],'Delimiter',',','preserveVariableNames',true);
% load divergent strain names (these are copied onto HD1, the rest onto HD2)
load('strainsList/divergent.mat','strains');
% initialise onHD1 variable to say if the file is on HD1
onHD1 = false(size(featureTable,1),1);
% go through the divergent strains to set their onHD1 status to true
for strainCtr = 1:numel(strains)
    strain = strains{strainCtr};
    strainLogInd = strcmp(featureTable.strain_name,strain);
    onHD1(strainLogInd) = true;
end
% append onHD1 variable to featureTable
featureTable.onHD1 = onHD1;

%% Go through each file that is present on this HD to extract features
% if HD1
%     fileInd = find(featureTable.onHD1);
% else
%     fileInd = find(~featureTable.onHD1);
% end
fileInd = 1:numel(featureTable.onHD1);

% initialise
filenames = repmat({''},numel(fileInd),1);
featVals = NaN(numel(fileInd),numel(featureBaseNames)*n_expandedFeats);
featNames = {};

for fileCtr = 1:numel(fileInd)
    disp(['Extracting features for File ' num2str(fileCtr) ' out of ' num2str(numel(fileInd)) ])
    %% Load data
    % get featuresN filename from Results folder
    filename = [strrep(strrep(char(featureTable.dirname(fileInd(fileCtr))),'/Volumes/behavgenom_archive$/Serena/','/Volumes/Ashur DT2/'),'MaskedVideos','Results'),...
        '/', strrep(char(featureTable.filename(fileInd(fileCtr))),'.hdf5','_featuresN.hdf5')];
    if ~featureTable.onHD1(fileInd(fileCtr))
        filename = strrep(filename,'DT2/','DT2 1/');
    end
    % load tracking data from featuresN.hdf5
    trajData = h5read(filename,'/trajectories_data');
    blobFeats = h5read(filename,'/blob_features');
    features = h5read(filename,'/timeseries_data');
    
    %% Subdivide and filter data
    % get worm category indices
    [singleWormLogInd,multiWormLogInd,clusterLogInd,pausedMwLogInd] = findClusters(trajData,blobFeats,clusterArea,features);
    % drop blobs that do not persist for more than 1 second
    tempBlobLogInd = findTempBlobs(trajData,blobFeats,frameRate);
    singleWormLogInd = singleWormLogInd & ~tempBlobLogInd;
    multiWormLogInd = multiWormLogInd & ~tempBlobLogInd;
    clusterLogInd = clusterLogInd & ~tempBlobLogInd;
    pausedMwLogInd = pausedMwLogInd & ~tempBlobLogInd;
    % get movie segment indices
    [phase1LogInd, phase2LogInd, phase3LogInd, notPhase1LogInd] = restrictPhase(trajData,frameRate);
    
    %% Calculate base features to expand and add to the main featuresTable
    [blobSpeed, d_blobSpeed, d_blobSpeed_abs] = calculateBlobSpeed(trajData, blobFeats,frameRate); % units in microns per second
    features2add = {blobSpeed,d_blobSpeed,d_blobSpeed_abs,blobFeats.area,blobFeats.compactness,blobFeats.perimeter,blobFeats.solidity,blobFeats.quirkiness,blobFeats.hu0,blobFeats.hu1,blobFeats.hu2,blobFeats.hu3,blobFeats.hu4,blobFeats.hu5,blobFeats.hu6,blobFeats.box_length,blobFeats.box_width,blobFeats.box_orientation};
    assert(numel(features2add) == numel(featureBaseNames),'The number of calculated base features does not match the number of feature base names.');
    
    %% Expand features
    for featureCtr = 1:numel(featureBaseNames)
        % Get base feature values and name
        feature = features2add{featureCtr};
        featureBaseName = featureBaseNames{featureCtr};
        % Expand features based on stats, worm type, and movie phase
        [values, names] = expandFeature(feature,featureBaseName,singleWormLogInd,multiWormLogInd,clusterLogInd,pausedMwLogInd,phase1LogInd,phase2LogInd,phase3LogInd,notPhase1LogInd);
        % Add results
        featVals(fileCtr,((featureCtr-1)*n_expandedFeats+1):featureCtr*n_expandedFeats) = values;
        if fileCtr==1 % only need to get featNames once
            featNames = horzcat(featNames,names);
        end
    end
    
    %% Get filename to join with the main features table later
    filenames{fileCtr} = char(featureTable.filename(fileInd(fileCtr))); % create this to join with the master FeaturesTable later
 end

%% Turn extracted features into a table
blobFeatureTable = array2table(featVals);
blobFeatureTable.Properties.VariableNames = featNames;
blobFeatureTable.filename = filenames;

%% Join new features table with the original
featureTable = outerjoin(featureTable, blobFeatureTable,'MergeKeys',true);

%% Save new featureTable with today's timestamp to append to existing featureTable later
timestamp=datetime('today','Format','yyyyMMdd');
writetable(featureTable,['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fortyWorm/fortyWormFeaturesTable_' extractStamp '_new_' char(timestamp) '.csv'])
save(['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fortyWorm/fortyWormFeaturesTable_' extractStamp '_new_' char(timestamp) '.mat'],'blobFeatureTable')