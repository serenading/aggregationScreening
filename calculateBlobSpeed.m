function [blobSpeeds, d_blobSpeeds, d_blobSpeeds_abs] = calculateBlobSpeed(trajData, blobFeats,frameRate,speedSmoothFactor,dT)

%% function calculates blobSpeed (smoothed over 1 second unless otherwise specified) 
% and d_blobSpeed (over dT window of 1/3 second unless otherwise specified, calculated from smoothed speeds)

%% Inputs:
% frameRate: scalar, usually = 25 for Phenix, fps
% speedSmoothFactor: scalar, the time window in seconds to smooth speed calculation over
% dT: scalar, the time window in seconds to calculate acceleration over

%% Outputs:
% blobSpeeds: unsigned speed (microns/s)
% d_blobSpeeds: signed acceleration (microns/s^2)
% d_blobSpeeds_abs: absolute value of acceleration (microns/s^2)

% check number of inputs and use default values if necessary
if nargin<5
    dT = 1/3; 
    if nargin<4
        speedSmoothFactor = 1;
        if nargin<3
            frameRate = 25;
        end
    end
end

% preallocate matrix to hold speed values
blobSpeeds = NaN(size(blobFeats.coord_x));
d_blobSpeeds = NaN(size(blobFeats.coord_x));
n_frames_speedsmooth = frameRate*speedSmoothFactor; % number of frames to smooth speed over
n_frames_dT = floor(dT*frameRate); % number of frames to calculate d_speed over

% go through each blob
uniqueBlobs = unique(trajData.worm_index_joined);
for blobCtr = 1:numel(uniqueBlobs)
    thisBlobLogInd = trajData.worm_index_joined == uniqueBlobs(blobCtr);
    % check that the blob exists for the smooth duration
    if nnz(thisBlobLogInd)>n_frames_speedsmooth
        blob_x = blobFeats.coord_x(thisBlobLogInd);
        blob_y = blobFeats.coord_y(thisBlobLogInd);
        blobSpeed = NaN(size(blob_x));
        % calculate blob speed
        blob_xBefore = vertcat(NaN(n_frames_speedsmooth,1),blob_x);
        blob_xAfter = vertcat(blob_x,NaN(n_frames_speedsmooth,1));
        dx = blob_xAfter-blob_xBefore;
        blob_yBefore = vertcat(NaN(n_frames_speedsmooth,1),blob_y);
        blob_yAfter = vertcat(blob_y,NaN(n_frames_speedsmooth,1));
        dy = blob_yAfter-blob_yBefore;
        blobSpeed = sqrt(dx.^2+dy.^2);
        blobSpeed = blobSpeed((n_frames_speedsmooth+1):end);
        % calculate blob acceleration
        speedsBefore = vertcat(NaN(n_frames_dT,1),blobSpeed);
        speedsAfter = vertcat(blobSpeed,NaN(n_frames_dT,1));
        d_blobSpeed = speedsAfter-speedsBefore;
        d_blobSpeed = d_blobSpeed((n_frames_dT+1):end);
        % check that indices for this blob is consecutive and that blobSpeed and d_blobSpeed are the same lengths
        assert(sum(diff(find(thisBlobLogInd))) == numel(diff(find(thisBlobLogInd))),'Indices for this blob are not consecutive.')
        assert(numel(d_blobSpeed) == numel(blobSpeed));
        % add blob speed to blobSpeeds matrix
        blobSpeeds(thisBlobLogInd) = blobSpeed;
        d_blobSpeeds(thisBlobLogInd) = d_blobSpeed;
    end
end

% convert blobSpeeds units into microns/second
blobSpeeds = blobSpeeds/speedSmoothFactor;
% convert d_blobSpeeds units into microns/second^2
d_blobSpeeds = d_blobSpeeds/dT;
% get absolute value of d_blobSpeeds
d_blobSpeeds_abs = abs(d_blobSpeeds);
end