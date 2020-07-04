clear
close all

clusterArea = 3;
frameRate = 25;

%% load
featureTablefilename = '/Users/sding/OneDrive - Imperial College London/aggScreening/results/fiveWorm/fiveWormFeaturesTable_20200519_153722.csv';
featureTable = readtable(featureTablefilename,'Delimiter',',','preserveVariableNames',true);

%% initialise and preallocate
% filenames
nfiles = numel(featureTable.filename);
filenames = repmat({''},nfiles,1);
% subdivision stats
divNames = {'tempBlobFraction','swFraction','mwFraction','clusterFraction','pausedMwFraction'};
divVals = NaN(nfiles,numel(divNames));

%% go through files
for fileCtr = 1:nfiles
    disp(['Extracting features for File ' num2str(fileCtr) ' out of ' num2str(nfiles) ])
    %% Load data
    % get featuresN filename from Results folder
    filename = [strrep(strrep(char(featureTable.dirname(fileCtr)),'/Volumes/behavgenom_archive$/Serena/','/Volumes/Ashur DT2/'),'MaskedVideos','Results'),...
        '/', strrep(char(featureTable.filename(fileCtr)),'.hdf5','_featuresN.hdf5')];
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
    
    % get div stats (can use "plotDivStats.m' to visualise this)
    nInd = numel(tempBlobLogInd);
    divVals(fileCtr,:) = [nnz(tempBlobLogInd),nnz(singleWormLogInd),nnz(multiWormLogInd),nnz(clusterLogInd),nnz(pausedMwLogInd)]/nInd; 

    %% Get filename to join with the main features table later
    filenames{fileCtr} = char(featureTable.filename(fileCtr)); % create this to join with the master FeaturesTable later
end

newFeatureTable = array2table(divVals);
newFeatureTable.Properties.VariableNames = divNames;
newFeatureTable.filename = filenames;

%% Save new featureTable with today's timestamp to append to existing featureTable later
timestamp=datetime('today','Format','yyyyMMdd');
save(['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fiveWorm/newFeatsToAdd/fiveWormFeaturesTable_20200519_153722_new_' char(timestamp) '.mat'],'newFeatureTable')