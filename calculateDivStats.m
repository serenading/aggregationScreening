clear
close all

addpath('auxiliary/')

%% Script extracts div statistics from tracking data to generate a feature table for joining onto the master feature table. Script is based on calculateBlobFeatures.m.
% author: serenading. July 2020.

%% Specify analysis parameters

% which HD is plugged in right now?
whichHD = 0; % Scalar. Specify 1 or 2. Specify 0 if both 1 and 2 are plugged in.

% which worm density?
wormNum = 40; % 40 or 5. 

% set which feature extraction timestamp to use.
extractStamp = '20200519_153722';

% number of expanded features per file, based on the expansion function.
n_expandedFeats = 119;

% some file attributes
frameRate = 25; % 25 frames per second, read by frameRate = double(h5readatt(strrep(filename,'featuresN','skeletons'),'/plate_worms','expected_fps'));

%% load feature table
[featureTable,fileInd] = extractHDLocation(wormNum,extractStamp,whichHD);

%% initialise and preallocate
filenames = repmat({''},numel(fileInd),1);
divStats = NaN(numel(fileInd),n_expandedFeats);

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
    
    %% Calculate and expand features
    [values, names] = expandDivFeature(singleWormLogInd,multiWormLogInd,clusterLogInd,pausedMwLogInd,tempBlobLogInd,onFoodLogInd,foodEdgeLogInd,offFoodLogInd,phase1LogInd,phase2LogInd,phase3LogInd,notPhase1LogInd);
   
    %% Write results
    divStats(fileCtr,:) = values;
    
    %% Get filename to join with the main features table later
    filenames{fileCtr} = char(featureTable.filename(fileInd(fileCtr))); % create this to join with the master FeaturesTable later
end

%% Turn extracted features into a table
newFeatureTable = array2table(divStats);
newFeatureTable.Properties.VariableNames = names;
newFeatureTable.filename = filenames;

%% Save new featureTable with today's timestamp to append to existing featureTable later
timestamp=datetime('today','Format','yyyyMMdd');
if wormNum ==40
    save(['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fortyWorm/newFeatsToAdd/fortyWormFeaturesTable_' extractStamp '_new_' char(timestamp) '_divstats.mat'],'newFeatureTable')
elseif wormNum ==5
    save(['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fiveWorm/newFeatsToAdd/fiveWormFeaturesTable_' extractStamp '_new_' char(timestamp) '_divstats.mat'],'newFeatureTable')
end
appendFeatsToFeatureTable(newFeatureTable,wormNum,extractStamp);