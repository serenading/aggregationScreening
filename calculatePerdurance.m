function [sw_perdurance,mw_perdurance,cluster_perdurance,sw_frameDist,mw_frameDist,cluster_frameDist] = calculatePerdurance(blobFeats,trajData,singleWormLogInd,multiWormLogInd,clusterLogInd,sw_frameDist,mw_frameDist,cluster_frameDist)

%% function calculates cluster perdurance (in number of frames)

sw_feature = blobFeats.perdurance(singleWormLogInd);
mw_feature = blobFeats.perdurance(multiWormLogInd);
cluster_feature = blobFeats.perdurance(clusterLogInd);

swUniqueInd = unique(sw_feature);
sw_perdurance = NaN(length(swUniqueInd),1);
for swCtr = 1:length(swUniqueInd)
    swInd = swUniqueInd(swCtr);
    sw_perdurance(swCtr) = numel(find(sw_feature == swInd));
end
swFrames = trajData.frame_number(singleWormLogInd)';
if ~isempty(swFrames)
    magicInd = diff([0 diff([swFrames]) 0]==1);
    continuousFrames = find(magicInd == -1) - find(magicInd == 1) + 1; % list lengths of consecutive frames
    singleFrames = length(swFrames) - sum(continuousFrames(:)); % find number of single frames
    sw_frameDist = [sw_frameDist ones(1,singleFrames) continuousFrames]; % compile distribution of consecutive frames by adding each movie
end

mwUniqueInd = unique(mw_feature);
mw_perdurance = NaN(length(mwUniqueInd),1);
for mwCtr = 1:length(mwUniqueInd)
    mwInd = mwUniqueInd(mwCtr);
    mw_perdurance(mwCtr) = numel(find(mw_feature == mwInd));
end
mwFrames = trajData.frame_number(multiWormLogInd)';
if ~isempty(mwFrames)
    magicInd = diff([0 diff([mwFrames]) 0]==1);
    continuousFrames = find(magicInd == -1) - find(magicInd == 1) + 1; % list lengths of consecutive frames
    singleFrames = length(mwFrames) - sum(continuousFrames(:)); % find number of single frames
    mw_frameDist = [mw_frameDist ones(1,singleFrames) continuousFrames]; % compile distribution of consecutive frames by adding each movie
end

clusterUniqueInd = unique(cluster_feature);
cluster_perdurance = NaN(length(clusterUniqueInd),1);
for clusterCtr = 1:length(clusterUniqueInd)
    clusterInd = clusterUniqueInd(clusterCtr);
    cluster_perdurance(clusterCtr) = numel(find(cluster_feature == clusterInd));
end
clusterFrames = trajData.frame_number(clusterLogInd)';
if ~isempty(clusterFrames)
    magicInd = diff([0 diff([clusterFrames]) 0]==1);
    continuousFrames = find(magicInd == -1) - find(magicInd == 1) + 1; % list lengths of consecutive frames
    singleFrames = length(clusterFrames) - sum(continuousFrames(:)); % find number of single frames
    cluster_frameDist = [cluster_frameDist ones(1,singleFrames) continuousFrames]; % compile distribution of consecutive frames by adding each movie
end