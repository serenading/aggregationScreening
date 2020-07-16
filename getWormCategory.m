function [singleWormLogInd,multiWormLogInd,clusterLogInd,pausedMwLogInd,tempBlobLogInd] = getWormCategory(trajData,blobFeats,tsFeatures,frameRate)

%% function generates clusterLogInd the same size as trajData files to detect clusters based on blob area.

%% INPUT: 
% trajData = h5read(filename,'/trajectories_data'); filename being the featuresN.hdf5.
% blobFeats = h5read(filename,'/blob_features');
% tsFeatures = h5read(filename,'/timeseries_data');
% frameRate = double(h5readatt(strrep(filename,'featuresN','skeletons'),'/plate_worms','expected_fps'));

%% OUTPUT:
% logical indices for single worms (blobs that skeletonise), multiworms
% (blobs that do not skeletonise), clusters (blobs that have min relative
% area), or paused multiworm (blobs that do not skeletonise and have zero
% motion_mode

% define the minimum size relative to a single worm area for a blob to be called a cluster; 
clusterArea = 3;

% find single and multi worms
singleWormLogInd = logical(trajData.was_skeletonized);
multiWormLogInd = ~logical(trajData.was_skeletonized);
% find clusters
swAreas = blobFeats.area(singleWormLogInd);
swArea = median(swAreas);
normAreas = blobFeats.area/swArea;
clusterLogInd =  normAreas > clusterArea;
% find paused multi worms
pausedMwLogInd = tsFeatures.motion_mode == 0 & trajData.was_skeletonized==0; % visual inspection of movies revealed this category of special interest

% drop blobs that do not persist for more than 1 second
tempBlobLogInd = findTempBlobs(trajData,blobFeats,frameRate);
singleWormLogInd = singleWormLogInd & ~tempBlobLogInd;
multiWormLogInd = multiWormLogInd & ~tempBlobLogInd;
clusterLogInd = clusterLogInd & ~tempBlobLogInd;
pausedMwLogInd = pausedMwLogInd & ~tempBlobLogInd;

end