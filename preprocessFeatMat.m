% Function pro-processes features matrix for PCA analysis.
% Input & Output: file-by-feature matrix containing feature values, class double.
% Output: column index of dropped features

function [featureMat,droppedCols] = preprocessFeatMat(featureMat)

% specify processing parameters
dropFeatThreshold = 0.2; % the maximum fraction of NaN values that a feature can have before being dropped
droppedCols = [];

% initial feature numbers
n_featsStart = size(featureMat,2);

% drop features with too many NaN's
featBefore = size(featureMat,2);
numNanFeat = sum(isnan(featureMat),1);
n_files = size(featureMat,1);
colsToDrop = numNanFeat > n_files*dropFeatThreshold; % get logical index for features with too many NaN's
featureMat = featureMat(:,~colsToDrop);
disp([ num2str(nnz(colsToDrop)) ' out of  ' num2str(featBefore) ' features dropped due to too many NaN values'])
droppedCols = [droppedCols find(colsToDrop)];

% impute nan values to global mean
featMeans = nanmean(featureMat);
imputeCtr = 0;
for featCtr = 1:size(featureMat, 2)
    nanInds = isnan(featureMat(:, featCtr));
    if nnz(nanInds)>0
        featureMat(nanInds, featCtr) = featMeans(featCtr);
        imputeCtr = imputeCtr+1;
    end
end
disp([ num2str(imputeCtr) ' out of  ' num2str(size(featureMat,2)) ' features have NaN values imputed'])

% drop features with zero standard deviation
featBefore = size(featureMat,2);
featStds = std(featureMat);
colsToDrop = featStds == 0 ; % get logical index for features with too many NaN's
featureMat = featureMat(:,~colsToDrop);
disp([ num2str(nnz(colsToDrop)) ' out of  ' num2str(featBefore) ' features dropped due to zero standard deviation'])
droppedCols = [droppedCols find(colsToDrop)];

% z-normalise feature matrix
featureMat = normalize(featureMat,1);
disp('z-normalisation applied to feature matrix')

% final feature numbers
n_featsEnd= size(featureMat,2);

% check feature number
assert(numel(droppedCols) == (n_featsStart - n_featsEnd));
disp(['Start feature number: ' num2str(n_featsStart) '; End feature number: ' num2str(n_featsEnd) '.'])
end