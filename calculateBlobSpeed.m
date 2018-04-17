function [blobSpeeds] = calculateBlobSpeed(trajData, blobFeats,frameRate,blobSpeedSmoothFactor)

%% function calculates blobSpeed (smoothed over 1 second unless otherwise specified) 
if nargin<4
    blobSpeedSmoothFactor = frameRate;
end

% preallocate matrix to hold speed values
blobSpeeds = NaN(size(blobFeats.coord_x));

% go through each blob
for blobCtr = 1:numel(unique(trajData.worm_index_joined))
    blobLogInd = trajData.worm_index_joined == blobCtr;
    % check that the blob lasts for at least 1 second
    if nnz(blobLogInd)>blobSpeedSmoothFactor
        % calculate blob speed
        blobInd = find(blobLogInd);
        blob_x = blobFeats.coord_x(blobInd);
        blob_y = blobFeats.coord_y(blobInd);
        blobSpeed = NaN(size(blob_x));
        for stepCtr = 1:(length(blob_x)-blobSpeedSmoothFactor)
            dx = blob_x(stepCtr+blobSpeedSmoothFactor)-blob_x(stepCtr);
            dy = blob_y(stepCtr+blobSpeedSmoothFactor)-blob_y(stepCtr);
            blobSpeed(stepCtr) = sqrt(dx.^2+dy.^2);
        end
        % add blob speed to blobSpeeds matrix
        blobSpeeds(blobInd(1):blobInd(end))=blobSpeed;
    end
end

% turn blobSpeeds units into microns/second
blobSpeeds = blobSpeeds/blobSpeedSmoothFactor*frameRate;

end