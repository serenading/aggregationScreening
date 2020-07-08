clear
close all

addpath('auxiliary/')

%% Script extracts Tierpsy blob features from 40 worm tracking data to generate a feature table for joining onto the master feature table

% author: serenading. June 2020.

%% Specify analysis parameters

% which HD is plugged in right now?
whichHD = 0; % Scalar. Specify 1 or 2. Specify 0 if both 1 and 2 are plugged in.

% which worm density?
wormNum = 5; % 40 or 5. 

% set which feature extraction timestamp to use (doesn't matter as we are calculating our own feats)
extractStamp = '20200519_153722';

% some file attributes
frameRate = 25; % 25 frames per second, read by frameRate = double(h5readatt(strrep(filename,'featuresN','skeletons'),'/plate_worms','expected_fps'));

% which base features to expand and add to features table? (this works for
% 40 worms only)
if wormNum ==40
    featureBaseNames = {'blobSpeed', 'd_blobSpeed', 'd_blobSpeed_abs','blobArea','blobPerimeter','blobCompactness','blobSolidity','blobQuirkiness','blobBoxLength','blobBoxWidth','blobOrientation','blobHu0','blobHu1','blobHu2','blobHu3','blobHu4','blobHu5','blobHu6'};
    n_expandedFeats = 80; % number of expanded features per base feature. = 100 for expandfeature, = 80 for expandfeature2
    % how many times the area of a single worm counts as a cluster?
    clusterArea = 3;
elseif wormNum ==5
    featureBaseNames = {};
    n_expandedFeats = NaN;
end

%% load feature table
if wormNum == 40
    %% find which files are on which HD
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
    if whichHD ==1
        fileInd = find(featureTable.onHD1);
    elseif whichHD ==2
        fileInd = find(~featureTable.onHD1);
    elseif whichHD ==0 % both HDs plugged in
        fileInd = 1:numel(featureTable.onHD1);
    end
elseif wormNum ==5
    %% all files are on HD1
    assert(whichHD == 1, ['HD' num2str(whichHD) ' is plugged in, but all five worm data are on HD1.']);
    featureTable = readtable(['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fiveWorm/fiveWormFeaturesTable_' extractStamp '.csv'],'Delimiter',',','preserveVariableNames',true);
    fileInd = 1:size(featureTable,1);
end

%% initialise and preallocate
% filenames
filenames = repmat({''},numel(fileInd),1);
% subdivision stats
if wormNum ==40
    divNames = {'tempBlobFraction','swFraction','mwFraction','clusterFraction','pausedMwFraction','food_region_inside_fraction','food_region_edge_fraction','food_region_outside_fraction'};
    % extracted and expanded stats
    featVals = NaN(numel(fileInd),numel(featureBaseNames)*n_expandedFeats);
    featNames = {};
elseif wormNum ==5
    divNames = {'tempBlobFraction','swFraction','mwFraction','clusterFraction','pausedMwFraction'};
end
divVals = NaN(numel(fileInd),numel(divNames));

for fileCtr = 1:numel(fileInd)
    disp(['Extracting features for File ' num2str(fileCtr) ' out of ' num2str(numel(fileInd)) ])
    %% Load data
    % get featuresN filename from Results folder
    filename = [strrep(strrep(char(featureTable.dirname(fileInd(fileCtr))),'/Volumes/behavgenom_archive$/Serena/','/Volumes/Ashur DT2/'),'MaskedVideos','Results'),...
        '/', strrep(char(featureTable.filename(fileInd(fileCtr))),'.hdf5','_featuresN.hdf5')];
    % replace HD path only if both HD are plugged in and 1 is plugged in before 2.
    if whichHD == 0 && ~featureTable.onHD1(fileInd(fileCtr))
        filename = strrep(filename,'DT2/','DT2 1/');
    end
    % load tracking data from featuresN.hdf5
    trajData = h5read(filename,'/trajectories_data');
    blobFeats = h5read(filename,'/blob_features');
    tsFeatures = h5read(filename,'/timeseries_data');
    
    %% Subdivide and filter data
    % get worm category indices
    [singleWormLogInd,multiWormLogInd,clusterLogInd,pausedMwLogInd] = findClusters(trajData,blobFeats,clusterArea,tsFeatures);
    % drop blobs that do not persist for more than 1 second
    tempBlobLogInd = findTempBlobs(trajData,blobFeats,frameRate);
    singleWormLogInd = singleWormLogInd & ~tempBlobLogInd;
    multiWormLogInd = multiWormLogInd & ~tempBlobLogInd;
    clusterLogInd = clusterLogInd & ~tempBlobLogInd;
    pausedMwLogInd = pausedMwLogInd & ~tempBlobLogInd;
    if wormNum ==40
        % get movie segment indices
        [phase1LogInd, phase2LogInd, phase3LogInd, notPhase1LogInd] = restrictPhase(trajData,frameRate);
        % get food region indices
        [onFoodLogInd,foodEdgeLogInd,offFoodLogInd] = getFoodRegion(tsFeatures);
        % get div stats (can use "plotDivStats.m' to visualise this)
        nInd = numel(tempBlobLogInd);
        divVals(fileCtr,:) = [nnz(tempBlobLogInd),nnz(singleWormLogInd),nnz(multiWormLogInd),nnz(clusterLogInd),nnz(pausedMwLogInd),...
            nnz(onFoodLogInd),nnz(foodEdgeLogInd),nnz(offFoodLogInd)]/nInd;
        
        %% Calculate base features to expand and add to the main featuresTable
        [blobSpeed, d_blobSpeed, d_blobSpeed_abs] = calculateBlobSpeed(trajData, blobFeats,frameRate); % units in microns per second
        features2add = {blobSpeed, d_blobSpeed, d_blobSpeed_abs,blobFeats.area,blobFeats.perimeter,blobFeats.compactness,blobFeats.solidity,blobFeats.quirkiness,blobFeats.box_length,blobFeats.box_width,blobFeats.box_orientation,blobFeats.hu0,blobFeats.hu1,blobFeats.hu2,blobFeats.hu3,blobFeats.hu4,blobFeats.hu5,blobFeats.hu6};
        assert(numel(features2add) == numel(featureBaseNames),'The number of calculated base features does not match the number of feature base names.');
        
        %% Expand features
        for featureCtr = 1:numel(featureBaseNames)
            % Get base feature values and name
            feature = features2add{featureCtr};
            featureBaseName = featureBaseNames{featureCtr};
            % Expand features based on stats, worm type, food region, and movie phase
            [values, names] = expandFeature2(feature,featureBaseName,singleWormLogInd,multiWormLogInd,clusterLogInd,pausedMwLogInd,onFoodLogInd,foodEdgeLogInd,offFoodLogInd);
            % Add results
            featVals(fileCtr,((featureCtr-1)*n_expandedFeats+1):featureCtr*n_expandedFeats) = values;
            if fileCtr==1 % only need to get featNames once
                featNames = horzcat(featNames,names);
            end
        end
    end
    %% Get filename to join with the main features table later
    filenames{fileCtr} = char(featureTable.filename(fileInd(fileCtr))); % create this to join with the master FeaturesTable later
end

%% Turn extracted features into a table
if wormNum==40
    newFeatureTable = array2table([divVals,featVals]);
    newFeatureTable.Properties.VariableNames = [divNames,featNames];
elseif wormNum==5
    newFeatureTable = array2table(divVals);
    newFeatureTable.Properties.VariableNames = divNames;
end
newFeatureTable.filename = filenames;

%% Save new featureTable with today's timestamp to append to existing featureTable later
timestamp=datetime('today','Format','yyyyMMdd');
if wormNum ==40
    save(['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fortyWorm/newFeatsToAdd/fortyWormFeaturesTable_' extractStamp '_new_' char(timestamp) '.mat'],'newFeatureTable')
elseif wormNum ==5
    save(['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fiveWorm/newFeatsToAdd/fiveWormFeaturesTable_' extractStamp '_new_' char(timestamp) '.mat'],'newFeatureTable')
end
appendFeatsToFeatureTable(newFeatureTable,wormNum,extractStamp);