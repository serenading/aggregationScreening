clear
close all

addpath('auxiliary/')

%% Script extracts Tierpsy blob features from 40 worm tracking data to generate a feature table for joining onto the master feature table

% author: serenading. June 2020.

% TODO: expand d_blob features; calculate and expand normalised blob features. 
% TODO: replace all d_blobSpeeds values. 
% TODO: automate feature expansion and complete feature expansion. combine
% the two expansion scripts now. 

%% Specify analysis parameters

% which HD is plugged in right now?
whichHD = 0; % Scalar. Specify 1 or 2. Specify 0 if both 1 and 2 are plugged in.

% which worm density?
wormNum = 40; % 40 or 5.

% set which feature extraction timestamp to use (doesn't matter as we are calculating our own feats)
extractStamp = '20200519_153722';

% some file attributes
frameRate = 25; % 25 frames per second, read by frameRate = double(h5readatt(strrep(filename,'featuresN','skeletons'),'/plate_worms','expected_fps'));

% which base features to expand and add to features table? (this works for
% 40 worms only)
featureBaseNames = {'blobSpeed', 'd_blobSpeed', 'd_blobSpeed_abs','blobArea','blobPerimeter','blobCompactness','blobSolidity','blobQuirkiness','blobBoxLength','blobBoxWidth','blobOrientation','blobHu0','blobHu1','blobHu2','blobHu3','blobHu4','blobHu5','blobHu6'};
n_expandedFeats = 80; % number of expanded features per base feature. = 100 for expandfeature, = 80 for expandfeature2

%% load feature table
[featureTable,fileInd] = extractHDLocation(40,extractStamp,whichHD);

%% initialise and preallocate
% filenames
filenames = repmat({''},numel(fileInd),1);
% extracted and expanded stats
featVals = NaN(numel(fileInd),numel(featureBaseNames)*n_expandedFeats);
featNames = {};

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
    [singleWormLogInd,multiWormLogInd,clusterLogInd,pausedMwLogInd,tempBlobLogInd] = getWormCategory(trajData,blobFeats,tsFeatures,frameRate);
    % get movie segment indices
    [phase1LogInd, phase2LogInd, phase3LogInd, notPhase1LogInd] = getMoviePhase(trajData,frameRate);
    % get food region indices
    [onFoodLogInd,foodEdgeLogInd,offFoodLogInd] = getFoodRegion(tsFeatures);
    
    %% Calculate base features to expand and add to the main featuresTable
    [blobSpeeds, d_blobSpeeds, d_blobSpeeds_abs, blobSpeeds_norm, d_blobSpeeds_norm, d_blobSpeeds_norm_abs] = ...
        calculateBlobSpeed(trajData,blobFeats,tsFeatures,frameRate); % units in microns per second
    features2add = {blobSpeed, d_blobSpeed, d_blobSpeed_abs, blobSpeeds_norm, d_blobSpeeds_norm, d_blobSpeeds_norm_abs,...
        blobFeats.area,blobFeats.perimeter,blobFeats.compactness,blobFeats.solidity,blobFeats.quirkiness,...
        blobFeats.box_length,blobFeats.box_width,blobFeats.box_orientation,...
        blobFeats.hu0,blobFeats.hu1,blobFeats.hu2,blobFeats.hu3,blobFeats.hu4,blobFeats.hu5,blobFeats.hu6};
    assert(numel(features2add) == numel(featureBaseNames),'The number of calculated base features does not match the number of feature base names.');
    
    %% Expand features
    for featureCtr = 1:numel(featureBaseNames)
        % Get base feature values and name
        feature = features2add{featureCtr};
        featureBaseName = featureBaseNames{featureCtr};
        % Expand features based on stats, worm type, food region, and movie phase
        [values, names] = expandBlobFeature2(feature,featureBaseName,singleWormLogInd,multiWormLogInd,clusterLogInd,pausedMwLogInd,onFoodLogInd,foodEdgeLogInd,offFoodLogInd);
        % Add results
        featVals(fileCtr,((featureCtr-1)*n_expandedFeats+1):featureCtr*n_expandedFeats) = values;
        if fileCtr==1 % only need to get featNames once
            featNames = horzcat(featNames,names);
        end
    end
    
    %% Expand with normalised features
    [values, names] = expandNormalisedBlobFeature;
    
    %% Get filename to join with the main features table later
    filenames{fileCtr} = char(featureTable.filename(fileInd(fileCtr))); % create this to join with the master FeaturesTable later
end

%% Turn extracted features into a table
newFeatureTable = array2table(featVals);
newFeatureTable.Properties.VariableNames = featNames;
newFeatureTable.filename = filenames;

%% Save new featureTable with today's timestamp to append to existing featureTable later
timestamp=datetime('today','Format','yyyyMMdd');
save(['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fortyWorm/newFeatsToAdd/fortyWormFeaturesTable_' extractStamp '_new_' char(timestamp) '_blobFeats.mat'],'newFeatureTable')
appendFeatsToFeatureTable(newFeatureTable,wormNum,extractStamp);

%%%%%%%%%%%%%%%%%%%%%%
%% LOCAL FUNCTION
%%%%%%%%%%%%%%%%%%%%%%

function [blobSpeeds, d_blobSpeeds, d_blobSpeeds_abs, blobSpeeds_norm, d_blobSpeeds_norm, d_blobSpeeds_norm_abs] = calculateBlobSpeed(trajData, blobFeats,tsFeatures,frameRate,speedSmoothFactor,dT)

%% function calculates blobSpeed (smoothed over 1 second unless otherwise specified) 
% and d_blobSpeed (over dT window of 1/3 second unless otherwise specified, calculated from smoothed speeds)

%% Inputs:
% frameRate: scalar, usually = 25 for Phenix, fps
% speedSmoothFactor: scalar, the time window in seconds to smooth speed calculation over
% dT: scalar, the time window in seconds to calculate acceleration over

%% Outputs:
% blobSpeeds: unsigned speed (microns/s)
% d_blobSpeeds: signed acceleration (microns/s^2)
% d_blobSpeeds_abs: absolute value of acceleration (microns/s^2)
% blobSpeeds_norm: unsigned speed normalised by worm lengths (microns/s/mm)
% d_blobSpeeds_norm: signed acceleration normalised by worm lengths
% d_blobSpeeds_norm_abs: absolute value of normalised acceleration

% check number of inputs and use default values if necessary
if nargin<5
    dT = 1/3; 
    if nargin<4
        speedSmoothFactor = 1;
        if nargin<3
            frameRate = 25;
        end
    end
end

% preallocate matrix to hold speed values
blobSpeeds = NaN(size(blobFeats.coord_x));
d_blobSpeeds = NaN(size(blobFeats.coord_x));
blobSpeeds_norm = NaN(size(blobFeats.coord_x));
d_blobSpeeds_norm = NaN(size(blobFeats.coord_x));
n_frames_speedsmooth = frameRate*speedSmoothFactor; % number of frames to smooth speed over
n_frames_dT = floor(dT*frameRate); % number of frames to calculate d_speed over

% go through each blob
uniqueBlobs = unique(trajData.worm_index_joined);
for blobCtr = 1:numel(uniqueBlobs)
    thisBlobLogInd = trajData.worm_index_joined == uniqueBlobs(blobCtr);
    % check that the blob exists for the smooth duration
    if nnz(thisBlobLogInd)>n_frames_speedsmooth
        blob_x = blobFeats.coord_x(thisBlobLogInd);
        blob_y = blobFeats.coord_y(thisBlobLogInd);
        blobSpeed = NaN(size(blob_x));
        % calculate blob speed
        blob_xBefore = vertcat(NaN(n_frames_speedsmooth,1),blob_x);
        blob_xAfter = vertcat(blob_x,NaN(n_frames_speedsmooth,1));
        dx = blob_xAfter-blob_xBefore;
        blob_yBefore = vertcat(NaN(n_frames_speedsmooth,1),blob_y);
        blob_yAfter = vertcat(blob_y,NaN(n_frames_speedsmooth,1));
        dy = blob_yAfter-blob_yBefore;
        blobSpeed = sqrt(dx.^2+dy.^2);
        blobSpeed = blobSpeed((n_frames_speedsmooth+1):end);
        blobSpeed = blobSpeed/speedSmoothFactor; % convert blobSpeed units into microns/second
        % normalise blob speed with worm length
        lengths = tsFeatures.length(logical(trajData.was_skeletonized)); % in microns
        medianLength = nanmedian(lengths)/1000; % median length in mm
        blobSpeed_norm = blobSpeed/medianLength;
        % calculate blob acceleration
        speedsBefore = vertcat(NaN(n_frames_dT,1),blobSpeed);
        speedsAfter = vertcat(blobSpeed,NaN(n_frames_dT,1));
        d_blobSpeed = speedsAfter-speedsBefore;
        d_blobSpeed = d_blobSpeed((n_frames_dT+1):end);
        d_blobSpeed = d_blobSpeed/dT; % convert d_blobSpeed units into microns/second^2
        % calculate blob acceleration based on normalised speed
        speedsBefore_norm = speedsBefore/medianLength;
        speedsAfter_norm = speedsAfter/medianLength;
        d_blobSpeed_norm = speedsAfter_norm - speedsBefore_norm;
        d_blobSpeed_norm = d_blobSpeed_norm((n_frames_dT+1):end);
        d_blobSpeed_norm = d_blobSpeed_norm/dT;
        % check that indices for this blob is consecutive and that blobSpeed and d_blobSpeed are the same lengths
        % (this step is slow, and previous running of the script shows the indices are consecutive for all files concerned)
        % assert(sum(diff(find(thisBlobLogInd))) == numel(diff(find(thisBlobLogInd))),'Indices for this blob are not consecutive.')
        assert(numel(d_blobSpeed) == numel(blobSpeed));
        % add blob speed to blobSpeeds matrix
        blobSpeeds(thisBlobLogInd) = blobSpeed;
        d_blobSpeeds(thisBlobLogInd) = d_blobSpeed;
        blobSpeeds_norm = blobSpeed_norm;
        d_blobSpeeds_norm(thisBlobLogInd) = d_blobSpeed_norm;
    end
end
% get absolute value of d_blobSpeeds
d_blobSpeeds_abs = abs(d_blobSpeeds);
d_blobSpeeds_norm_abs = abs(d_blobSpeeds_norm);
end