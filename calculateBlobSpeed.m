function blobSpeeds = calculateBlobSpeed(trajData, blobFeats,frameRate,pixelToMicron,speedSmoothFactorInSec)

%% function calculates blobSpeed (smoothed over 1 second unless otherwise specified) 
if nargin<5
    speedSmoothFactorInSec = 1;
end

% preallocate matrix to hold speed values
blobSpeeds = NaN(size(blobFeats.coord_x));

% go through each blob
uniqueBlobs = unique(trajData.worm_index_joined);
for blobCtr = 1:numel(uniqueBlobs)
    thisBlobLogInd = trajData.worm_index_joined == uniqueBlobs(blobCtr);
    % check that the blob exists for the smooth duration
    if nnz(thisBlobLogInd)>speedSmoothFactorInSec*frameRate
        % calculate blob speed
        blob_x = blobFeats.coord_x(thisBlobLogInd);
        blob_y = blobFeats.coord_y(thisBlobLogInd);
        blobSpeed = NaN(size(blob_x));
        for stepCtr = 1:(numel(blob_x)-speedSmoothFactorInSec*frameRate)
            dx = blob_x(stepCtr+speedSmoothFactorInSec*frameRate)-blob_x(stepCtr);
            dy = blob_y(stepCtr+speedSmoothFactorInSec*frameRate)-blob_y(stepCtr);
            blobSpeed(stepCtr) = sqrt(dx.^2+dy.^2);
        end
        % check that indices for this blob is consecutive
        assert(sum(diff(find(thisBlobLogInd))) == numel(diff(find(thisBlobLogInd))),'Indices for this blob are not consecutive.')
        % add blob speed to blobSpeeds matrix
        blobSpeeds(thisBlobLogInd)=blobSpeed;
    end
end

% turn blobSpeeds units into microns/second
blobSpeeds = blobSpeeds/speedSmoothFactorInSec*pixelToMicron;

end