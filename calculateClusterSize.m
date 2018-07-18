function [mwSizeNorm,clusterSizeNorm] = calculateClusterSize(blobFeats,trajData,frameLogInd,numFrames2Sample,multiWormLogInd,clusterLogInd,medianSwFeat)

% phase restrict movies
validFrames = trajData.frame_number(frameLogInd);
% randomly sample a number of frames
sampleFrames = datasample(validFrames,numFrames2Sample,'Replace',false);
% preallocate variables
clusterSizeNorm = [];
mwSizeNorm = [];
% go through each frame
for frameCtr = 1:numel(sampleFrames)
    frame = sampleFrames(frameCtr);
    % grow variable with normalised cluster sizes from each frame
    mwInd = find(multiWormLogInd & trajData.frame_number == frame);
    if ~isempty(mwInd)
        mwSizeNorm = vertcat(mwSizeNorm,blobFeats.clusterSize(mwInd)/medianSwFeat);
    end
    clusterInd = find(clusterLogInd & trajData.frame_number == frame);
    if ~isempty(clusterInd)
        clusterSizeNorm = vertcat(clusterSizeNorm,blobFeats.clusterSize(clusterInd)/medianSwFeat);
    end
end