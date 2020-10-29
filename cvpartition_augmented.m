function holdoutCVP = cvpartition_augmented(classLabels)

%% Function divides augmented features up into training and test datasets, so that one entire replicate is held out so it is not at all seen in training.
% Input:
% classLabels: nx1 cell array containing classification response variables where n is the number of rows in the features table to be partitioned.
% Output:
% holdoutCVP: structure containing two nx1 logical vectors for test and training sets (training set is held out).
% Author: @serenading. Oct 2020.

augFactor = 5; % each experimental replicate is augmented to produce five rows of features. 

% Get a list of unique strains for classification
uniqueStrains = unique(classLabels);
% Preallocate
holdoutFirstInd = NaN(numel(uniqueStrains),1);
% Go through each strain
for strainCtr = 1:numel(uniqueStrains)
    % Get strain name
    strain = uniqueStrains(strainCtr);
    % Get the row index corresponding to this strain
    strainInd = find(strcmp(classLabels,strain));
    n_expReps = numel(strainInd)/augFactor;
    % Check that there are exactly 5x n_expReps number of strainInd
    assert(floor(n_expReps) == n_expReps)
    % Choose a random replicate out of the two or three for each of the strains
    holdoutRepIdx = randi([1,n_expReps]);
    % Record the first row index of the chosen replicate
    holdoutFirstInd(strainCtr) = strainInd(augFactor*(holdoutRepIdx-1)+1);
end
% Expand the holdout indices to include all five augmented rows per experimental replicate
holdoutInd = [holdoutFirstInd,holdoutFirstInd+1,holdoutFirstInd+2,holdoutFirstInd+3,holdoutFirstInd+4];
holdoutInd = sort(holdoutInd(:));
% Get logical index for holdout set
holdoutLogInd = false(numel(classLabels),1);
holdoutLogInd(holdoutInd) = true;
% Check that the holdout set encapsulates one experimental replicate and five augmented result rows per strain
assert(numel(unique(classLabels(holdoutLogInd))) == numel(uniqueStrains));
assert(nnz(holdoutLogInd) == numel(uniqueStrains)*augFactor);
% Write logical index into structure akin to cvpartition
holdoutCVP.training = ~holdoutLogInd;
holdoutCVP.test = holdoutLogInd;

end