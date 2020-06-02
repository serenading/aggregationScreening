%% Function pro-processes features matrix with several steps: 
% 1. Drop features with zero standard deviation.
% 2. Drop features with too many NaN's.
% 3. Impute nan values to global mean.
% 4. z-normalise feature matrix.

% Input & Output: file-by-feature matrix containing feature values, class double.

% Output: column index of dropped features

function [featureMat,dropLogInd] = preprocessFeatMat(featureMat)

% specify processing parameters
dropFeatThreshold = 0.2; % the maximum fraction of NaN values that a feature can have before being dropped

% initialise
n_featsStart = size(featureMat,2);
dropLogInd = false(1,n_featsStart); % this keeps track of all dropped features

% drop features with zero standard deviation
featBefore = size(featureMat,2);
featStds = nanstd(featureMat);
colsToDrop = featStds == 0 ; % get logical index for features with zero standard deviation
dropLogInd(colsToDrop) = true;
disp([ num2str(nnz(colsToDrop)) ' out of  ' num2str(featBefore) ' features dropped due to zero standard deviation.'])

% drop features with too many NaN's
featBefore = size(featureMat,2);
numNanFeat = sum(isnan(featureMat),1);
n_files = size(featureMat,1);
colsToDrop = numNanFeat > n_files*dropFeatThreshold; % get logical index for features with too many NaN's
dropLogInd(colsToDrop) = true;
disp([ num2str(nnz(colsToDrop)) ' out of  ' num2str(featBefore) ' features dropped due to too many NaN values'...
    ' (more than ' num2str(dropFeatThreshold*100) '% NaN).'])
featureMat = featureMat(:,~dropLogInd); % drop the feats 

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
disp([ num2str(imputeCtr) ' out of  ' num2str(size(featureMat,2)) ' features have NaN values imputed.'])

% z-normalise feature matrix
featureMat = normalize(featureMat,1);
disp('z-normalisation applied to feature matrix.')

% final feature numbers
n_featsEnd= size(featureMat,2);

% check feature number
assert(nnz(dropLogInd) == (n_featsStart - n_featsEnd));
disp(['Feature matrix pre-processing: start feature number: ' num2str(n_featsStart) '; end feature number: ' num2str(n_featsEnd) '.'])

end