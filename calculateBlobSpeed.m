function [blobSpeeds] = calculateBlobSpeed(trajData, blobFeats,frameRate,pixelToMicron,speedSmoothFactorInSec)

%% function calculates blobSpeed (smoothed over 1 second unless otherwise specified) 
if nargin<5
    speedSmoothFactorInSec = 1;
end

% preallocate matrix to hold speed values
blobSpeeds = NaN(size(blobFeats.coord_x));

% go through each blob
for blobCtr = 1:numel(unique(trajData.worm_index_joined))
    blobLogInd = trajData.worm_index_joined == blobCtr;
    % check that the blob exists for the smooth duration
    if nnz(blobLogInd)>speedSmoothFactorInSec*frameRate
        % calculate blob speed
        blobInd = find(blobLogInd);
        blob_x = blobFeats.coord_x(blobInd);
        blob_y = blobFeats.coord_y(blobInd);
        blobSpeed = NaN(size(blob_x));
        for stepCtr = 1:(length(blob_x)-speedSmoothFactorInSec*frameRate)
            dx = blob_x(stepCtr+speedSmoothFactorInSec*frameRate)-blob_x(stepCtr);
            dy = blob_y(stepCtr+speedSmoothFactorInSec*frameRate)-blob_y(stepCtr);
            blobSpeed(stepCtr) = sqrt(dx.^2+dy.^2);
        end
        % add blob speed to blobSpeeds matrix
        blobSpeeds(blobInd(1):blobInd(end))=blobSpeed;
    end
end

% turn blobSpeeds units into microns/second
blobSpeeds = blobSpeeds/speedSmoothFactorInSec*pixelToMicron;

end