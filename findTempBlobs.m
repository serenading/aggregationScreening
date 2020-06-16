function tempBlobLogInd = findTempBlobs(trajData,blobFeats,frameRate)

%% function generates tempBlobLogInd the same size as trajData files to detect blobs that do not persist for more than 1 second.

% initialise logical index
tempBlobLogInd = false(size(trajData.worm_index_joined));
% unique blobs 
uniqueBlobs = unique(trajData.worm_index_joined);
% go through each blob
for blobCtr = 1:numel(uniqueBlobs)
    thisBlobLogInd = trajData.worm_index_joined == uniqueBlobs(blobCtr);
    % check that the blob persists for at least one second
    if nnz(thisBlobLogInd)<=frameRate
        % find the indices of this blob
        thisBlobInd = find(thisBlobLogInd);
        % turn those indices false
        tempBlobLogInd(thisBlobInd)=true;
    end
end