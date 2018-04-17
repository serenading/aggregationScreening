function [tempBlobLogInd] = findTempBlobs(trajData,blobFeats,frameRate)

%% function generates tempBlobLogInd the same size as trajData files to detect blobs that do not persist for more than 1 second.

% initialise logical index
tempBlobLogInd = true(size(trajData.worm_index_joined));
% go through each blob
for blobCtr = 1:numel(unique(trajData.worm_index_joined))
    blobLogInd = trajData.worm_index_joined == blobCtr;
    % check that the blob persists for at least second
    if nnz(blobLogInd)>=frameRate
        % find the indices of this blob
        blobInd = find(blobLogInd);
        % turn those indices false
        tempBlobLogInd(blobInd(1):blobInd(end))=false;
    end
end