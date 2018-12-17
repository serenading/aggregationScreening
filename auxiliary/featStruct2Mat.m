function [featMat, wormInds, featNames] = featStruct2Mat(fileName, minLength, featFlag)

% FEATSTRUCT2MAT imports the mean features from a tierpsy HDF5 feature file
% and converts them into a feature matrix.  Features derived from
% trajectories less than minLength are not included.
% 
% Input
%   fileName  - the full path of the file to be imported.
%   minLength - the minimum trajectory length (in frames) required for data
%               to be included.
%   featFlag  - logical vector of length 726 indicating whether each mean 
%               feature is to be included. Mean feature table has 730
%               entries, but first 4 are metadata not features.
% 
% Output
%   featMat   - a numTrajectories x numFeatures matrix of mean feature
%               values. May contain NaNs for features that could not be
%               calculated.
%   wormInds  - the indices used to identify each trajectory. Can be used
%               to get skeleton and trajectory data corresponding to
%               mean feature data.
%   featNames - the names of each of the features

% check inputs
if ~islogical(featFlag)
    error('featFlag must be a vector of logicals.')
end

% load feature data
try
    featData = h5read(fileName, '/features_means');
catch
    warning('Input file cannot be read. Possibly empty.')
    featMat = [];
    wormInds = [];
    featNames = [];
    return
end

% get worm indices
wormInds = featData.worm_index;

% get trajectory lengths (in frames)
trajLengths = featData.n_frames;

% get feature names
featNames = fieldnames(featData);
featNames(1:4) = []; % drop metadata fields

% check featFlag length
if size(featNames, 1) ~= length(featFlag)
    error('featFlag must have the same length as number of features.')
end

% drop unused feature names
featNames = featNames(featFlag);

% convert to cell, dropping metadata entries at start
featCell = struct2cell(featData);
featCell(1:4) = [];

% loop through features to reshape into numTrajectories x numFeatures
% matrix
keepInds = find(featFlag);
featMat = NaN(size(featCell{1}, 1), numel(keepInds));
for ii = 1:numel(keepInds)
    featMat(:, ii) = featCell{keepInds(ii)};
end

% drop short trajectories
featMat = featMat(trajLengths > minLength, :);