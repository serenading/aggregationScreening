function [singleWormLogInd,multiWormLogInd,clusterLogInd,pausedMwLogInd] = findClusters(trajData,blobFeats,clusterArea,features)

%% function generates clusterLogInd the same size as trajData files to detect clusters based on blob area.

%% INPUT: 
% trajData = h5read(filename,'/trajectories_data'); filename being the featuresN.hdf5.
% blobFeats = h5read(filename,'/blob_features');
% clusterArea: minimum cluster relative area compared to a single worm area;
% features = h5read(filename,'/timeseries_data');

%% OUTPUT:
% logical indices for single worms (blobs that skeletonise), multiworms
% (blobs that do not skeletonise), clusters (blobs that have min relative
% area), or paused multiworm (blobs that do not skeletonise and have zero
% motion_mode


% find single and multi worms
singleWormLogInd = logical(trajData.was_skeletonized);
multiWormLogInd = ~logical(trajData.was_skeletonized);
% find clusters
swAreas = blobFeats.area(singleWormLogInd);
swArea = median(swAreas);
normAreas = blobFeats.area/swArea;
clusterLogInd =  normAreas > clusterArea;
% find paused multi worms
pausedMwLogInd = features.motion_mode == 0 & trajData.was_skeletonized==0; % visual inspection of movies revealed this category of special interest