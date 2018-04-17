function [singleWormLogInd,multiWormLogInd,clusterLogInd] = findClusters(trajData,blobFeats,clusterArea)

%% function generates clusterLogInd the same size as trajData files to detect clusters based on blob area.

%% INPUT: 
%% trajData = h5read(filename,'/trajectories_data');
%% blobFeats = h5read(filename,'/blob_features');
%% clusterArea: minimum cluster relative area compared to a single worm area. 4 by default if unspecified.

%% OUTPUT:
%% logical indices for single worms (blobs that skeletonise), multiworms (blobs that do not skeletonise), and clusters (blobs that have min relative area of 3.5 or as specified)

if nargin<3
    clusterArea = 4; % if unspecified, then use 4 as the default clusterArea
end
    

singleWormLogInd = logical(trajData.is_good_skel);
multiWormLogInd = ~logical(trajData.is_good_skel);
swAreas = blobFeats.area(singleWormLogInd);
swArea = median(swAreas);
normAreas = blobFeats.area/swArea;
clusterLogInd =  normAreas > clusterArea;
